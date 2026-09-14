# 12_direct_ml.R  (V2 — fixes look-ahead via feature lagging)
#
# Direct Parametric Portfolio Policy (PPP) following Brandt, Santa-Clara,
# and Valkanov (2009). End-to-end policy: features at time t-1 -> equity
# weight at time t, no intermediate regime classification.
#
# CRITICAL FIX vs V1:
#   V1 used X(t) to predict w_eq(t), then R_p(t) = w_eq(t)*r_eq(t). Since
#   X(t) contained `eq_ret` (= r_eq(t)), the optimizer learned a perfect-
#   timing rule with look-ahead bias. Result: RR_median=1.215, unrealistic.
#
# V2 timing convention:
#   w_eq(t) = 0.7 * sigmoid(X(t-1)' * beta)
#   R_p(t)  = w_eq(t) * r_eq(t) + (1-w_eq(t)) * r_bond(t)
#
# This matches Framing 1 (regime-based ML): features observed at t-1 decide
# allocation that earns return realized at t. On DNB scenarios, this means
# year-t features become input for year-(t+1) decision.
#
# Implementation:
#   - Historical training: lag ALL features by 1 month (was: contemporaneous)
#   - DNB inference: features at year t produce w_eq for year t+1
#   - For year 1 of DNB, use warm-start from Dec 2025 historical features
#
# Other choices unchanged from V1:
#   - Layer 1 linear policy: w_eq = 0.7 * sigmoid(X' * beta)
#   - Objective C3 (Tu-Zhou utility): mean - (gamma/2) * var - lambda * ||beta||^2
#   - gamma = 5, eq_cap = 0.70 (IORP II)
#   - 7-fold time-series CV (W_init=84, V_win=24, gap=12)
#   - L-BFGS-B optimization with analytical gradient
#
# Input/Output: same as V1.

library(dplyr)
library(tidyr)
library(readxl)

set.seed(42)

# ---- 1. Configuration ---------------------------------------------------
GAMMA       <- 5
EQ_CAP      <- 0.7
N_FOLDS     <- 7
W_INIT      <- 84
V_WIN       <- 24
GAP         <- 12
LAMBDA_GRID <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50)

# ---- 2. Load inputs ------------------------------------------------------
hist         <- readRDS("Data/features_historical.rds")
feature_cols <- readRDS("Data/feature_cols.rds")

allocations  <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")

y1  <- as.matrix(readRDS("Data/yield_1y.rds"))  * 100
y5  <- as.matrix(readRDS("Data/yield_5y.rds"))  * 100
y10 <- as.matrix(readRDS("Data/yield_10y.rds")) * 100
y20 <- as.matrix(readRDS("Data/yield_20y.rds")) * 100
hicp_yoy_dnb <- as.matrix(read_excel(
  "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx",
  sheet = "5_EU_Price_Inflation", col_names = FALSE)) * 100

metrics_old <- readRDS("Data/metrics_full.rds")
wp_old      <- readRDS("Data/wealth_paths.rds")
mdd_old     <- readRDS("Data/max_drawdowns.rds")
to_old      <- readRDS("Data/turnovers.rds")
RR_old      <- readRDS("Data/replacement_ratios_with_hybrid_naive.rds")

# If V1 was previously appended, strip Direct_PPP rows to avoid duplication
if ("Direct_PPP" %in% metrics_old$strategy) {
  cat("Note: removing previous Direct_PPP rows from metrics/wealth/RR/MDD/turnover\n\n")
  metrics_old <- metrics_old %>% filter(strategy != "Direct_PPP")
  wp_old[["Direct_PPP"]] <- NULL
  mdd_old <- mdd_old[, colnames(mdd_old) != "Direct_PPP", drop = FALSE]
  to_old  <- to_old [, colnames(to_old)  != "Direct_PPP", drop = FALSE]
  RR_old  <- RR_old [, colnames(RR_old)  != "Direct_PPP", drop = FALSE]
}

N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

# ---- 3. Build historical bond returns ------------------------------------
hist <- hist %>%
  arrange(date) %>%
  mutate(
    bond_ret_m = (lag(y10) / 100) / 12 - 9 * (y10 - lag(y10)) / 100
  )

# ---- 4. CRITICAL FIX: lag ALL features by 1 month -----------------------
# X(t-1) predicts w_eq(t); R_p(t) uses r_eq(t), r_bond(t).
# Some features in feature_cols are ALREADY lagged (eq_ret_l1, etc); we lag
# them again by 1, making them effectively l2/l4/l13. This is acceptable:
# the point is to have a feature set entirely observable BEFORE the
# allocation decision.

hist_lagged <- hist %>%
  mutate(across(all_of(feature_cols), ~ lag(., 1), .names = "{.col}_lag1"))

# Now keep only rows where: lagged features are all non-NA AND bond return non-NA
lag_cols <- paste0(feature_cols, "_lag1")
hist_use <- hist_lagged %>%
  filter(if_all(all_of(lag_cols), ~ !is.na(.))) %>%
  filter(!is.na(bond_ret_m))

n_hist <- nrow(hist_use)
cat("Historical sample after lagging features:", n_hist, "monthly obs\n")
cat("  Date range:", as.character(min(hist_use$date)), "to",
    as.character(max(hist_use$date)), "\n\n")

# Standardize lagged features
X_full <- as.matrix(hist_use[, lag_cols])
colnames(X_full) <- feature_cols   # rename back to original for downstream use
feat_mean <- colMeans(X_full)
feat_sd   <- apply(X_full, 2, sd)
X_std <- scale(X_full, center = feat_mean, scale = feat_sd)

# Realised returns (contemporaneous with the decision they fund)
r_eq   <- hist_use$eq_ret
r_bond <- hist_use$bond_ret_m

# Design matrix with intercept
X_design <- cbind(intercept = 1, X_std)
p <- ncol(X_design)
cat("Design matrix:", n_hist, "obs x", p, "cols (incl. intercept)\n\n")

# ---- 5. Sanity check on what the top feature would have been ---------
# Without lagging, V1's top beta was eq_ret with +1.13 (current return =
# best predictor of current return). With lagging, eq_ret_lag1 is the
# previous-month return; its predictive power for next-month return is
# weak (momentum coefficient typically <0.05 in monthly returns).
cat("Univariate correlations of lagged features with r_eq (sanity):\n")
cors <- sapply(seq_len(ncol(X_std)),
               function(j) cor(X_std[, j], r_eq))
names(cors) <- feature_cols
cors_sorted <- sort(abs(cors), decreasing = TRUE)
for (j in 1:5) {
  nm <- names(cors_sorted)[j]
  cat(sprintf("  %-20s  |cor| = %.3f  (sign: %+.3f)\n",
              nm, cors_sorted[j], cors[nm]))
}
cat("\n(If top |cor| > 0.5 we still have look-ahead. Expected: all < 0.20.)\n\n")

# ---- 6. Objective + analytical gradient ---------------------------------
sigmoid <- function(z) {
  ifelse(z >= 0, 1 / (1 + exp(-z)), exp(z) / (1 + exp(z)))
}

portfolio_returns <- function(w_eq, r_eq, r_bond) {
  w_eq * r_eq + (1 - w_eq) * r_bond
}

neg_utility <- function(beta, X, r_eq, r_bond, gamma, eq_cap, lambda) {
  z <- as.numeric(X %*% beta)
  w_eq <- eq_cap * sigmoid(z)
  r_p <- portfolio_returns(w_eq, r_eq, r_bond)
  util <- mean(r_p) - (gamma / 2) * var(r_p)
  reg  <- lambda * sum(beta^2)
  -(util - reg)
}

neg_utility_grad <- function(beta, X, r_eq, r_bond, gamma, eq_cap, lambda) {
  z <- as.numeric(X %*% beta)
  sig <- sigmoid(z)
  dsig_dz <- sig * (1 - sig)
  w_eq <- eq_cap * sig
  spread <- r_eq - r_bond
  dRp_dbeta <- eq_cap * dsig_dz * spread
  dRp_dbeta_mat <- X * dRp_dbeta
  Tn <- length(r_eq)
  r_p <- portfolio_returns(w_eq, r_eq, r_bond)
  r_p_dev <- r_p - mean(r_p)
  grad_mean <- colMeans(dRp_dbeta_mat)
  grad_var  <- (2 / Tn) * colSums(r_p_dev * dRp_dbeta_mat)
  grad_reg  <- 2 * lambda * beta
  -(grad_mean - (gamma / 2) * grad_var - grad_reg)
}

# ---- 7. Time-series folds (same as Framing 1) ---------------------------
build_folds <- function(n, W_init, V_win, gap, K) {
  folds <- list()
  i <- 1
  start <- W_init + gap + 1
  while (length(folds) < K) {
    end <- start + V_win - 1
    if (end > n) break
    folds[[i]] <- list(
      train_idx = 1:(start - gap - 1),
      val_idx   = start:end
    )
    start <- start + V_win
    i <- i + 1
  }
  folds
}

folds <- build_folds(n_hist, W_INIT, V_WIN, GAP, N_FOLDS)
cat("Built", length(folds), "CV folds\n\n")

# ---- 8. Fit PPP ----------------------------------------------------------
fit_ppp <- function(X_train, r_eq_train, r_bond_train, lambda) {
  init_beta <- rep(0, ncol(X_train))
  res <- tryCatch(
    optim(
      par     = init_beta,
      fn      = neg_utility,
      gr      = neg_utility_grad,
      X       = X_train,
      r_eq    = r_eq_train,
      r_bond  = r_bond_train,
      gamma   = GAMMA,
      eq_cap  = EQ_CAP,
      lambda  = lambda,
      method  = "BFGS",
      control = list(maxit = 200, reltol = 1e-9)
    ),
    error = function(e) NULL
  )
  if (is.null(res) || res$convergence != 0) return(NULL)
  res$par
}

# ---- 9. Cross-validation ------------------------------------------------
cat("Running 7-fold CV across", length(LAMBDA_GRID), "lambda values...\n")
cv_results <- expand.grid(
  fold = seq_along(folds), lambda = LAMBDA_GRID,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
cv_results$val_util <- NA_real_

for (i in seq_len(nrow(cv_results))) {
  f      <- cv_results$fold[i]
  lambda <- cv_results$lambda[i]
  tr <- folds[[f]]$train_idx
  vl <- folds[[f]]$val_idx
  beta_hat <- fit_ppp(X_design[tr, ], r_eq[tr], r_bond[tr], lambda)
  if (is.null(beta_hat)) next
  z_val   <- as.numeric(X_design[vl, ] %*% beta_hat)
  w_val   <- EQ_CAP * sigmoid(z_val)
  rp_val  <- portfolio_returns(w_val, r_eq[vl], r_bond[vl])
  cv_results$val_util[i] <- mean(rp_val) - (GAMMA / 2) * var(rp_val)
}

lambda_summary <- cv_results %>%
  group_by(lambda) %>%
  summarise(mean_util = mean(val_util, na.rm = TRUE),
            n_valid = sum(!is.na(val_util)), .groups = "drop")
cat("\nCV utility by lambda:\n")
print(lambda_summary)

best_lambda <- lambda_summary$lambda[which.max(lambda_summary$mean_util)]
cat(sprintf("\nBest lambda = %.4f\n\n", best_lambda))

# Diagnostic: is best_lambda at an interior point (good) or at corners (bad)?
if (best_lambda == min(LAMBDA_GRID)) {
  cat("WARNING: best_lambda is at the LOWER bound. Possible underfitting of regularization.\n\n")
}
if (best_lambda == max(LAMBDA_GRID)) {
  cat("WARNING: best_lambda is at the UPPER bound. Possible no-signal scenario.\n\n")
}

# ---- 10. Final fit on full historical sample ----------------------------
cat("Fitting final PPP on full historical sample...\n")
beta_final <- fit_ppp(X_design, r_eq, r_bond, best_lambda)
if (is.null(beta_final)) stop("Final PPP fit failed.")

beta_named <- setNames(beta_final, colnames(X_design))
cat("\nFinal beta coefficients (top 10 by absolute value):\n")
top_abs <- names(sort(abs(beta_named), decreasing = TRUE))[1:10]
for (nm in top_abs) {
  cat(sprintf("  %-20s  beta = %+.4f\n", nm, beta_named[nm]))
}

# In-sample equity weight stats
z_in <- as.numeric(X_design %*% beta_final)
w_in <- EQ_CAP * sigmoid(z_in)
cat(sprintf("\nIn-sample equity weight statistics:\n"))
cat(sprintf("  min:    %.3f\n",   min(w_in)))
cat(sprintf("  max:    %.3f\n",   max(w_in)))
cat(sprintf("  mean:   %.3f\n",   mean(w_in)))
cat(sprintf("  median: %.3f\n",   median(w_in)))
cat(sprintf("  sd:     %.3f\n\n", sd(w_in)))

if (mean(w_in > EQ_CAP - 0.005) > 0.5) {
  cat("WARNING: w_eq binds the cap in >50% of months (degenerate policy).\n\n")
}
if (mean(w_in < 0.05) > 0.5) {
  cat("WARNING: w_eq near 0 in >50% of months (degenerate policy).\n\n")
}

# ---- 11. OOF predictions on historical -------------------------------
cat("Computing OOF predictions on historical sample...\n")
oof_w <- rep(NA_real_, n_hist)
for (f in seq_along(folds)) {
  tr <- folds[[f]]$train_idx
  vl <- folds[[f]]$val_idx
  beta_hat <- fit_ppp(X_design[tr, ], r_eq[tr], r_bond[tr], best_lambda)
  if (is.null(beta_hat)) next
  z_vl <- as.numeric(X_design[vl, ] %*% beta_hat)
  oof_w[vl] <- EQ_CAP * sigmoid(z_vl)
}
n_oof <- sum(!is.na(oof_w))
cat("  OOF coverage:", n_oof, "of", n_hist, "\n\n")

# ---- 12. Build DNB features and apply PPP ------------------------------
# CRITICAL: features for year t use info AT END OF YEAR t-1 (lag-1
# convention). For year 1, we warm-start from Dec 2025 historical features
# (the last row of hist).

hist_last <- tail(hist, 1)

build_features_year_lagged <- function(t_career) {
  # Returns features observed at END OF YEAR (t_career-1).
  # When t_career == 1: features come from Dec 2025 (historical).
  # When t_career >= 2: features come from DNB year (t_career-1).
  if (t_career == 1) {
    fv <- as.numeric(hist_last[, feature_cols])
    X_t <- matrix(fv, nrow = N_SCEN, ncol = length(feature_cols), byrow = TRUE)
    colnames(X_t) <- feature_cols
    return(X_t)
  }
  
  tp <- t_career - 1
  y1_t  <- y1[,  tp + 1]
  y5_t  <- y5[,  tp + 1]
  y10_t <- y10[, tp + 1]
  y20_t <- y20[, tp + 1]
  if (tp == 1) {
    y1_prev  <- y1[, 1]
    y10_prev <- y10[, 1]
  } else {
    y1_prev  <- y1[,  tp]
    y10_prev <- y10[, tp]
  }
  dy1_t  <- 100 * (y1_t  - y1_prev)
  dy10_t <- 100 * (y10_t - y10_prev)
  
  eq_ret_t   <- eq_ret[, tp]
  hicp_yoy_t <- hicp_yoy_dnb[, tp]
  
  slope_10_1_t   <- y10_t - y1_t
  slope_20_5_t   <- y20_t - y5_t
  curvature_t    <- 2 * y5_t - y1_t - y10_t
  real_yield10_t <- y10_t - hicp_yoy_t
  eq_drawdown_12m_t <- pmin(0, eq_ret_t)
  
  if (tp == 1) {
    eq_ret_l1_t    <- rep(hist_last$eq_ret,       N_SCEN)
    eq_ret_l3_t    <- rep(hist_last$eq_ret_l3,    N_SCEN)
    eq_ret_l12_t   <- rep(hist_last$eq_ret_l12,   N_SCEN)
    hicp_yoy_l3_t  <- rep(hist_last$hicp_yoy_l3,  N_SCEN)
    hicp_yoy_l12_t <- rep(hist_last$hicp_yoy_l12, N_SCEN)
    slope_l3_t     <- rep(hist_last$slope_l3,     N_SCEN)
    slope_l12_t    <- rep(hist_last$slope_l12,    N_SCEN)
    y10_l12_t      <- rep(hist_last$y10_l12,      N_SCEN)
  } else {
    eq_ret_l1_t   <- eq_ret[, tp - 1]
    eq_ret_l3_t   <- eq_ret_l1_t
    eq_ret_l12_t  <- eq_ret[, tp - 1]
    hicp_yoy_l3_t <- hicp_yoy_t
    hicp_yoy_l12_t <- hicp_yoy_dnb[, tp - 1]
    y1p <- y1[, tp]; y10p <- y10[, tp]
    sp  <- y10p - y1p
    slope_l3_t  <- sp
    slope_l12_t <- sp
    y10_l12_t   <- y10p
  }
  
  slope_x_infl_t <- slope_10_1_t * hicp_yoy_t
  realy_x_eq_t   <- real_yield10_t * eq_ret_l1_t
  
  cbind(eq_ret = eq_ret_t, hicp_yoy = hicp_yoy_t,
        y1 = y1_t, y5 = y5_t, y10 = y10_t, y20 = y20_t,
        dy1 = dy1_t, dy10 = dy10_t,
        slope_10_1 = slope_10_1_t, slope_20_5 = slope_20_5_t,
        curvature = curvature_t, real_yield10 = real_yield10_t,
        eq_drawdown_12m = eq_drawdown_12m_t,
        eq_ret_l1 = eq_ret_l1_t, eq_ret_l3 = eq_ret_l3_t,
        eq_ret_l12 = eq_ret_l12_t,
        hicp_yoy_l3 = hicp_yoy_l3_t, hicp_yoy_l12 = hicp_yoy_l12_t,
        slope_l3 = slope_l3_t, slope_l12 = slope_l12_t,
        y10_l12 = y10_l12_t,
        slope_x_infl = slope_x_infl_t, realy_x_eq = realy_x_eq_t)
}

cat("Applying PPP to 20,000 DNB scenarios with lagged features...\n")
w_direct <- matrix(NA_real_, nrow = N_SCEN, ncol = N_YEARS)
for (t in 1:N_YEARS) {
  X_t_raw <- build_features_year_lagged(t)
  X_t_std <- scale(X_t_raw, center = feat_mean, scale = feat_sd)
  X_t_design <- cbind(intercept = 1, X_t_std)
  z_t <- as.numeric(X_t_design %*% beta_final)
  w_direct[, t] <- EQ_CAP * sigmoid(z_t)
}

cat(sprintf("DNB equity weight statistics:\n"))
cat(sprintf("  min:    %.3f\n",   min(w_direct)))
cat(sprintf("  max:    %.3f\n",   max(w_direct)))
cat(sprintf("  mean:   %.3f\n",   mean(w_direct)))
cat(sprintf("  median: %.3f\n",   median(w_direct)))
cat(sprintf("  sd:     %.3f\n\n", sd(w_direct)))

cat("Year-by-year mean equity weight (first 10 years):\n")
print(round(colMeans(w_direct[, 1:10]), 3))
cat("\n")

# ---- 13. Wealth, MDD, turnover, RR -----------------------------------
compute_wealth_path <- function(w_eq) {
  wp <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wp[, t + 1] <- (wp[, t] + contrib) * (1 + port_ret)
  }
  wp
}

compute_mdd <- function(wp_matrix) {
  N <- nrow(wp_matrix)
  out <- numeric(N)
  for (i in 1:N) {
    path <- wp_matrix[i, ]
    rm <- cummax(path)
    dd <- ifelse(rm > 0, (path - rm) / rm, 0)
    out[i] <- min(dd)
  }
  out
}

compute_turnover <- function(w_eq) {
  Tn <- ncol(w_eq)
  if (Tn < 2) return(rep(0, nrow(w_eq)))
  rowSums(abs(w_eq[, 2:Tn] - w_eq[, 1:(Tn - 1)]))
}

cat("Computing wealth, MDD, turnover, RR for Direct_PPP...\n")
wp_direct  <- compute_wealth_path(w_direct)
mdd_direct <- compute_mdd(wp_direct)
to_direct  <- compute_turnover(w_direct)
RR_direct  <- wp_direct[, N_YEARS + 1] / af / salary_at_retirement

cat(sprintf("  RR_median = %.3f\n",   median(RR_direct)))
cat(sprintf("  RR_mean   = %.3f\n",   mean(RR_direct)))
cat(sprintf("  RR_sd     = %.3f\n",   sd(RR_direct)))
cat(sprintf("  Shortfall(70%%)  = %.3f\n", mean(RR_direct < 0.70)))
cat(sprintf("  CVaR_5%%          = %.3f\n",
            mean(RR_direct[RR_direct <= quantile(RR_direct, 0.05)])))
cat(sprintf("  MDD_mean        = %.3f\n", mean(mdd_direct)))
cat(sprintf("  Turnover_mean   = %.3f\n\n", mean(to_direct)))

# ---- 14. Append to metrics structures --------------------------------
metrics_row <- data.frame(
  strategy         = "Direct_PPP",
  RR_mean          = round(mean(RR_direct), 3),
  RR_median        = round(median(RR_direct), 3),
  RR_sd            = round(sd(RR_direct), 3),
  RR_p10           = round(quantile(RR_direct, 0.10), 3),
  RR_p25           = round(quantile(RR_direct, 0.25), 3),
  RR_p75           = round(quantile(RR_direct, 0.75), 3),
  RR_p90           = round(quantile(RR_direct, 0.90), 3),
  RR_iqr           = round(quantile(RR_direct, 0.75) - quantile(RR_direct, 0.25), 3),
  P_shortfall_50pct = round(mean(RR_direct < 0.50), 3),
  P_shortfall_60pct = round(mean(RR_direct < 0.60), 3),
  P_shortfall_70pct = round(mean(RR_direct < 0.70), 3),
  P_shortfall_80pct = round(mean(RR_direct < 0.80), 3),
  CVaR_5pct        = round(mean(RR_direct[RR_direct <= quantile(RR_direct, 0.05)]), 3),
  CVaR_10pct       = round(mean(RR_direct[RR_direct <= quantile(RR_direct, 0.10)]), 3),
  RR_sharpe_like   = round(mean(RR_direct) / sd(RR_direct), 3),
  MDD_mean         = round(mean(mdd_direct), 3),
  MDD_median       = round(median(mdd_direct), 3),
  MDD_p10          = round(quantile(mdd_direct, 0.10), 3),
  turnover_mean    = round(mean(to_direct), 3),
  turnover_median  = round(median(to_direct), 3),
  stringsAsFactors = FALSE,
  row.names = NULL
)

metrics_new <- bind_rows(metrics_old, metrics_row)

wp_new <- wp_old
wp_new[["Direct_PPP"]] <- wp_direct

mdd_new <- cbind(mdd_old, Direct_PPP = mdd_direct)
to_new  <- cbind(to_old,  Direct_PPP = to_direct)
RR_new  <- cbind(RR_old,  Direct_PPP = RR_direct)

# ---- 15. Comparison display ---------------------------------------
cat("\n--- METRICS COMPARISON: Direct_PPP vs baselines ---\n")
key_strategies <- c("60/40", "1/N", "Glide", "MV_plain", "RF",
                    "HybridNaive", "Direct_PPP")
print(metrics_new %>%
        filter(strategy %in% key_strategies) %>%
        select(strategy, RR_median, RR_sd, P_shortfall_70pct,
               CVaR_5pct, MDD_mean, RR_sharpe_like, turnover_mean),
      row.names = FALSE)

# ---- 16. Save all ------------------------------------------------
saveRDS(list(
  beta = beta_final, beta_named = beta_named,
  lambda = best_lambda, lambda_grid_results = lambda_summary,
  feat_mean = feat_mean, feat_sd = feat_sd,
  feature_names = colnames(X_design),
  oof_w = oof_w, folds = folds,
  version = "V2_lagged_features"
), "Data/direct_ml_model.rds")

saveRDS(w_direct, "Data/direct_ml_allocations.rds")
saveRDS(metrics_new, "Data/metrics_full.rds")
saveRDS(wp_new,      "Data/wealth_paths.rds")
saveRDS(mdd_new,     "Data/max_drawdowns.rds")
saveRDS(to_new,      "Data/turnovers.rds")
saveRDS(RR_new,      "Data/replacement_ratios_with_hybrid_naive.rds")

cat("\nSaved (V2, look-ahead corrected):\n")
cat("  Data/direct_ml_model.rds\n")
cat("  Data/direct_ml_allocations.rds\n")
cat("  Data/metrics_full.rds          (10 strategies)\n")
cat("  Data/wealth_paths.rds          (10 strategies)\n")
cat("  Data/max_drawdowns.rds         (10 strategies)\n")
cat("  Data/turnovers.rds             (10 strategies)\n")
cat("  Data/replacement_ratios_with_hybrid_naive.rds (10 strategies)\n")
