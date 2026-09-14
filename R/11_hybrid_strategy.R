# 11_hybrid_strategy.R  (V2 — OOF ML predictions, eliminates look-ahead)
#
# Input:  Data/features_historical.rds       (for historical returns calibration)
#         Data/features_with_regimes.rds     (for historical regime labels)
#         Data/ml_cv_results.rds             (OOF predictions from CV)  ← NEW
#         Data/dnb_equity_allocations.rds    (for ML allocations on DNB)
#         Data/equity_returns.rds            (DNB)
#         Data/bond_returns_10y.rds          (DNB)
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/mv_weights.rds
#         Data/replacement_ratios_all.rds    (8 baseline strategies)
#
# Output: Data/replacement_ratios_with_hybrid.rds  (≥10 strategies)
#         Data/hybrid_weights.rds                  (combination weights)
#         Data/robustness_hybrid_representatives.rds
#
# V2 CHANGE: ML historical returns use OUT-OF-FOLD predictions from
# 04_train_ml.R cross-validation, NOT perfect regime labels. This
# eliminates the in-sample look-ahead bias of V1, producing a genuine
# Sharpe calibration.
#
# Methodological note (Option X): the calibration sample is the 168 monthly
# observations that have OOF predictions for all 3 ML models (Lasso, Ridge,
# RF). Rule-based and MV historical returns are also restricted to the same
# 168 obs for comparability. Evaluation occurs on 20,000 DNB scenarios,
# which are out-of-sample with respect to both training and weight
# calibration.
#
# Reference: Tu & Zhou (2011) — combination logic applied to strategy-
# family diversification (model-free, model-based, machine-learning).

library(dplyr)
library(tidyr)

# ---- 1. Load inputs ------------------------------------------------------
hist          <- readRDS("Data/features_historical.rds")
allocations   <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret        <- readRDS("Data/equity_returns.rds")
bd_ret        <- readRDS("Data/bond_returns_10y.rds")
salary        <- readRDS("Data/salary_path.rds")
af            <- readRDS("Data/annuity_factor.rds")
mv_weights    <- readRDS("Data/mv_weights.rds")
RR_baseline   <- readRDS("Data/replacement_ratios_all.rds")
cv_res        <- readRDS("Data/ml_cv_results.rds")  # NEW: OOF predictions

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

# ---- 2. Reconstruct equity-weight matrices for all 8 base strategies ----
w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
w_glide <- matrix(rep(seq(0.80, 0.30, length.out = N_YEARS),
                      each = N_SCEN), N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_weights$plain["equity"]), N_SCEN, N_YEARS)
w_mv_lw <- matrix(as.numeric(mv_weights$lw["equity"]),    N_SCEN, N_YEARS)
w_lasso <- allocations[, , "Lasso"]
w_ridge <- allocations[, , "Ridge"]
w_rf    <- allocations[, , "RF"]

# ---- 3. Wealth and RR helpers --------------------------------------------
compute_wealth_final <- function(w_eq) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

compute_rr <- function(w_eq) {
  wf <- compute_wealth_final(w_eq)
  wf / af / salary_at_retirement
}

summarize_rr <- function(RR_vec, label) {
  data.frame(
    strategy        = label,
    RR_mean         = round(mean(RR_vec), 3),
    RR_median       = round(median(RR_vec), 3),
    RR_sd           = round(sd(RR_vec), 3),
    P_shortfall_70  = round(mean(RR_vec < 0.70), 3),
    CVaR_5pct       = round(mean(RR_vec[RR_vec <= quantile(RR_vec, 0.05)]), 3),
    stringsAsFactors = FALSE
  )
}

# ---- 4. Build historical strategy returns (monthly, OOF-restricted) -----
# Bond returns proxy
hist <- hist %>%
  mutate(bond_ret_m = (lag(y10) / 100) / 12 - 9 * (y10 - lag(y10)) / 100)

n_full <- nrow(hist)
cat("Full historical sample:", n_full, "monthly obs\n")

# OOF mask: TRUE where Lasso, Ridge AND RF all have non-NA predictions
oof_mask <- !is.na(cv_res$oof_class$lasso) &
  !is.na(cv_res$oof_class$ridge) &
  !is.na(cv_res$oof_class$rf)
cat("OOF-valid obs (all 3 models predict):", sum(oof_mask), "\n")

# The hist data should align with cv_res; verify by checking dims
stopifnot(length(oof_mask) == n_full)

# Subset hist to OOF-valid rows
hist_oof <- hist[oof_mask, ]
hist_oof <- hist_oof %>% filter(!is.na(bond_ret_m))
oof_mask_aligned <- which(oof_mask & !is.na(hist$bond_ret_m))
n_oof <- length(oof_mask_aligned)
cat("OOF-valid obs (with valid bond return):", n_oof, "\n")
cat("Date range:", as.character(min(hist_oof$date)), "to",
    as.character(max(hist_oof$date)), "\n\n")

# Allocation map
alloc_map <- c(Normal_expansion = 0.70, Inflationary_exp = 0.55,
               Late_cycle = 0.45, Stagflation_risk = 0.30)

# ---- 5. Equity weight vectors per strategy on OOF subset -----------------
# Rule-based: constant
hw_6040  <- rep(0.60, n_oof)
hw_1N    <- rep(0.50, n_oof)
hw_glide <- seq(0.80, 0.30, length.out = n_oof)
# MV
hw_mv_p  <- rep(as.numeric(mv_weights$plain["equity"]), n_oof)
hw_mv_lw <- rep(as.numeric(mv_weights$lw["equity"]), n_oof)
# ML: USE OOF predictions instead of perfect regime labels
hw_lasso <- alloc_map[cv_res$oof_class$lasso[oof_mask_aligned]]
hw_ridge <- alloc_map[cv_res$oof_class$ridge[oof_mask_aligned]]
hw_rf    <- alloc_map[cv_res$oof_class$rf[oof_mask_aligned]]

# Sanity check: ML weights should be different across the 3 models in some months
cat("ML OOF weight distribution:\n")
cat("  Lasso unique values:", sort(unique(hw_lasso)), "\n")
cat("  Ridge unique values:", sort(unique(hw_ridge)), "\n")
cat("  RF unique values:   ", sort(unique(hw_rf)),    "\n\n")

cat("ML OOF agreement rates:\n")
cat("  Lasso == Ridge:", round(mean(hw_lasso == hw_ridge), 3), "\n")
cat("  Lasso == RF:   ", round(mean(hw_lasso == hw_rf),    3), "\n")
cat("  Ridge == RF:   ", round(mean(hw_ridge == hw_rf),    3), "\n\n")

# Historical monthly portfolio returns
hist_port_ret <- list(
  `60/40`    = hw_6040  * hist_oof$eq_ret + (1 - hw_6040)  * hist_oof$bond_ret_m,
  `1/N`      = hw_1N    * hist_oof$eq_ret + (1 - hw_1N)    * hist_oof$bond_ret_m,
  Glide      = hw_glide * hist_oof$eq_ret + (1 - hw_glide) * hist_oof$bond_ret_m,
  MV_plain   = hw_mv_p  * hist_oof$eq_ret + (1 - hw_mv_p)  * hist_oof$bond_ret_m,
  MV_LW      = hw_mv_lw * hist_oof$eq_ret + (1 - hw_mv_lw) * hist_oof$bond_ret_m,
  Lasso      = hw_lasso * hist_oof$eq_ret + (1 - hw_lasso) * hist_oof$bond_ret_m,
  Ridge      = hw_ridge * hist_oof$eq_ret + (1 - hw_ridge) * hist_oof$bond_ret_m,
  RF         = hw_rf    * hist_oof$eq_ret + (1 - hw_rf)    * hist_oof$bond_ret_m
)

cat("Historical OOF Sharpe ratios (annualised):\n")
for (s in names(hist_port_ret)) {
  r <- hist_port_ret[[s]]
  sr <- mean(r) / sd(r) * sqrt(12)
  cat(sprintf("  %-10s  Sharpe = %.3f\n", s, sr))
}
cat("\n")

# ---- 6. Calibration: Sharpe-optimal weights ------------------------------
calibrate_sharpe <- function(r1, r2, r3, step = 0.05) {
  best_sharpe <- -Inf
  best_w <- c(NA, NA, NA)
  grid <- seq(0, 1, by = step)
  for (w1 in grid) {
    for (w2 in seq(0, 1 - w1, by = step)) {
      w3 <- 1 - w1 - w2
      if (w3 < 0 || w3 > 1) next
      r_combo <- w1 * r1 + w2 * r2 + w3 * r3
      sr <- mean(r_combo) / sd(r_combo) * sqrt(12)
      if (sr > best_sharpe) {
        best_sharpe <- sr
        best_w <- c(w1, w2, w3)
      }
    }
  }
  list(weights = best_w, sharpe = best_sharpe)
}

# ---- 7. Hybrid combinations to test --------------------------------------
combinations <- list(
  default = list(
    mf = "60/40", mb = "MV_plain", dd = "RF",
    label = "60/40 + MV_plain + RF"
  ),
  alt_1N = list(
    mf = "1/N", mb = "MV_plain", dd = "RF",
    label = "1/N + MV_plain + RF"
  ),
  alt_Glide = list(
    mf = "Glide", mb = "MV_plain", dd = "RF",
    label = "Glide + MV_plain + RF"
  ),
  alt_MVLW = list(
    mf = "60/40", mb = "MV_LW", dd = "RF",
    label = "60/40 + MV_LW + RF"
  ),
  alt_Lasso = list(
    mf = "60/40", mb = "MV_plain", dd = "Lasso",
    label = "60/40 + MV_plain + Lasso"
  ),
  alt_Ridge = list(
    mf = "60/40", mb = "MV_plain", dd = "Ridge",
    label = "60/40 + MV_plain + Ridge"
  )
)

dnb_weights <- list(
  `60/40` = w_6040, `1/N` = w_1N, Glide = w_glide,
  MV_plain = w_mv_p, MV_LW = w_mv_lw,
  Lasso = w_lasso, Ridge = w_ridge, RF = w_rf
)

# ---- 8. Compute Naive and Calibrated Hybrid -----------------------------
hybrid_results <- list()
hybrid_weights_log <- list()

for (cmb_name in names(combinations)) {
  cmb <- combinations[[cmb_name]]
  cat(sprintf("\n=== Combination: %s ===\n", cmb$label))
  
  r1 <- hist_port_ret[[cmb$mf]]
  r2 <- hist_port_ret[[cmb$mb]]
  r3 <- hist_port_ret[[cmb$dd]]
  
  # Naive
  w_naive <- c(1/3, 1/3, 1/3)
  r_naive <- w_naive[1] * r1 + w_naive[2] * r2 + w_naive[3] * r3
  sharpe_naive <- mean(r_naive) / sd(r_naive) * sqrt(12)
  
  # Calibrated
  calib <- calibrate_sharpe(r1, r2, r3)
  w_calib <- calib$weights
  
  cat(sprintf("  Naive weights:      (%.2f, %.2f, %.2f), Sharpe = %.3f\n",
              w_naive[1], w_naive[2], w_naive[3], sharpe_naive))
  cat(sprintf("  Calibrated weights: (%.2f, %.2f, %.2f), Sharpe = %.3f\n",
              w_calib[1], w_calib[2], w_calib[3], calib$sharpe))
  
  # Apply to DNB
  dnb_w1 <- dnb_weights[[cmb$mf]]
  dnb_w2 <- dnb_weights[[cmb$mb]]
  dnb_w3 <- dnb_weights[[cmb$dd]]
  
  w_dnb_naive <- w_naive[1] * dnb_w1 + w_naive[2] * dnb_w2 + w_naive[3] * dnb_w3
  w_dnb_calib <- w_calib[1] * dnb_w1 + w_calib[2] * dnb_w2 + w_calib[3] * dnb_w3
  
  RR_naive <- compute_rr(w_dnb_naive)
  RR_calib <- compute_rr(w_dnb_calib)
  
  label_naive <- sprintf("HybridNaive_%s", cmb_name)
  label_calib <- sprintf("HybridCalib_%s", cmb_name)
  
  hybrid_results[[label_naive]] <- summarize_rr(RR_naive, label_naive) %>%
    mutate(combination = cmb$label, type = "Naive (1/3)",
           hist_sharpe = round(sharpe_naive, 3))
  hybrid_results[[label_calib]] <- summarize_rr(RR_calib, label_calib) %>%
    mutate(combination = cmb$label, type = "Calibrated (Sharpe)",
           hist_sharpe = round(calib$sharpe, 3))
  
  hybrid_weights_log[[cmb_name]] <- list(
    label = cmb$label,
    components = c(cmb$mf, cmb$mb, cmb$dd),
    w_naive = w_naive,
    w_calib = w_calib,
    sharpe_naive = sharpe_naive,
    sharpe_calib = calib$sharpe
  )
}

# ---- 9. Master table -----------------------------------------------------
hybrid_master <- bind_rows(hybrid_results)
baseline_summary <- lapply(colnames(RR_baseline), function(s) {
  summarize_rr(RR_baseline[, s], s) %>%
    mutate(combination = NA, type = "Baseline", hist_sharpe = NA)
}) %>% bind_rows()
master <- bind_rows(baseline_summary, hybrid_master)

# ---- 10. Display key results --------------------------------------------
cat("\n\n=== HYBRID DEFAULT (60/40 + MV_plain + RF) ===\n")
default_rows <- master %>%
  filter(strategy %in% c("60/40", "MV_plain", "RF",
                         "HybridNaive_default", "HybridCalib_default")) %>%
  select(strategy, RR_median, RR_mean, RR_sd,
         P_shortfall_70, CVaR_5pct, hist_sharpe)
print(default_rows, row.names = FALSE)

cat("\n=== ALL STRATEGIES (DEFAULT HYBRID) ===\n")
key_cols <- master %>%
  filter(type %in% c("Baseline", "Naive (1/3)", "Calibrated (Sharpe)")
         | strategy %in% c("HybridNaive_default", "HybridCalib_default")) %>%
  select(strategy, type, RR_median, RR_mean,
         P_shortfall_70, CVaR_5pct, hist_sharpe)
print(key_cols, row.names = FALSE)

cat("\n=== HYBRID ROBUSTNESS ACROSS REPRESENTATIVES ===\n")
robust_table <- master %>%
  filter(type %in% c("Naive (1/3)", "Calibrated (Sharpe)")) %>%
  select(combination, type, RR_median, P_shortfall_70, CVaR_5pct, hist_sharpe) %>%
  arrange(combination, type)
print(robust_table, row.names = FALSE)

# ---- 11. Save ------------------------------------------------------------
saveRDS(master,                "Data/replacement_ratios_with_hybrid.rds")
saveRDS(hybrid_weights_log,    "Data/hybrid_weights.rds")
saveRDS(hybrid_master,         "Data/robustness_hybrid_representatives.rds")
cat("\nSaved Data/replacement_ratios_with_hybrid.rds\n")
cat("Saved Data/hybrid_weights.rds\n")
cat("Saved Data/robustness_hybrid_representatives.rds\n")

# 08_metrics_v2.R
#
# Extends the baseline metrics pipeline to include the Naive Hybrid strategy
# (60/40 + MV_plain + RF, equal 1/3 weights) as the 9th strategy.
#
# Input:  Data/metrics_full.rds          (8 strategies baseline)
#         Data/wealth_paths.rds          (8 strategies)
#         Data/max_drawdowns.rds         (8 strategies)
#         Data/turnovers.rds             (8 strategies)
#         Data/dnb_equity_allocations.rds
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/mv_weights.rds
#         Data/replacement_ratios_all.rds
#
# Output: Same files, overwritten, with 9th strategy appended.
#         Data/replacement_ratios_with_hybrid_naive.rds (RR matrix 20000 x 9)
#
# Note: HybridCalib is excluded from main reporting (rejected for
# methodological reasons - see 11_hybrid_strategy.R V2).

library(dplyr)
library(tidyr)

# ---- 1. Load existing outputs --------------------------------------------
metrics_old   <- readRDS("Data/metrics_full.rds")
wp_old        <- readRDS("Data/wealth_paths.rds")
mdd_old       <- readRDS("Data/max_drawdowns.rds")
to_old        <- readRDS("Data/turnovers.rds")
RR_old        <- readRDS("Data/replacement_ratios_all.rds")

allocations   <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret        <- readRDS("Data/equity_returns.rds")
bd_ret        <- readRDS("Data/bond_returns_10y.rds")
salary        <- readRDS("Data/salary_path.rds")
af            <- readRDS("Data/annuity_factor.rds")
mv_weights    <- readRDS("Data/mv_weights.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

cat("Existing strategies in metrics_full:\n")
print(metrics_old$strategy)
cat("Adding: HybridNaive_default\n\n")

# ---- 2. Build HybridNaive equity-weight matrix --------------------------
w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_weights$plain["equity"]), N_SCEN, N_YEARS)
w_rf    <- allocations[, , "RF"]

w_hybrid <- (1/3) * w_6040 + (1/3) * w_mv_p + (1/3) * w_rf

cat("HybridNaive equity weight statistics:\n")
cat(sprintf("  min:    %.3f\n", min(w_hybrid)))
cat(sprintf("  max:    %.3f\n", max(w_hybrid)))
cat(sprintf("  mean:   %.3f\n", mean(w_hybrid)))
cat(sprintf("  median: %.3f\n\n", median(w_hybrid)))

# ---- 3. Wealth path (matrix 20000 x 43) ---------------------------------
# Year 0 = 0 wealth, year t = wealth at end of year t after contribution
# and return. Matches the convention of wealth_paths.rds.

# First verify the shape of an existing wealth path
existing_wp_shape <- dim(wp_old[[1]])
cat("Existing wealth_paths element shape:", existing_wp_shape, "\n")
# Should be 20000 x (N_YEARS+1) = 20000 x 43

wealth_path_hybrid <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)
for (t in 1:N_YEARS) {
  contrib  <- CONTRIB_RATE * salary[, t]
  port_ret <- w_hybrid[, t] * eq_ret[, t] +
    (1 - w_hybrid[, t]) * bd_ret[, t]
  wealth_path_hybrid[, t + 1] <-
    (wealth_path_hybrid[, t] + contrib) * (1 + port_ret)
}

wealth_final_hybrid <- wealth_path_hybrid[, N_YEARS + 1]
RR_hybrid <- wealth_final_hybrid / af / salary_at_retirement

cat("HybridNaive RR statistics:\n")
cat(sprintf("  mean:   %.3f\n", mean(RR_hybrid)))
cat(sprintf("  median: %.3f\n", median(RR_hybrid)))
cat(sprintf("  sd:     %.3f\n\n", sd(RR_hybrid)))

# ---- 4. Maximum drawdown per scenario -----------------------------------
# MDD = max drawdown of wealth trajectory per scenario
compute_mdd_per_scenario <- function(wp_matrix) {
  N <- nrow(wp_matrix)
  mdd_vec <- numeric(N)
  for (i in 1:N) {
    path <- wp_matrix[i, ]
    running_max <- cummax(path)
    # Drawdown defined relative to running peak (avoid div by 0 at t=0)
    dd <- ifelse(running_max > 0, (path - running_max) / running_max, 0)
    mdd_vec[i] <- min(dd)
  }
  mdd_vec
}

cat("Computing MDD for HybridNaive...\n")
mdd_hybrid <- compute_mdd_per_scenario(wealth_path_hybrid)
cat(sprintf("  MDD mean:   %.3f\n", mean(mdd_hybrid)))
cat(sprintf("  MDD median: %.3f\n\n", median(mdd_hybrid)))

# ---- 5. Turnover per scenario -------------------------------------------
# Turnover = sum over t of |w_eq[i, t] - w_eq[i, t-1]|
compute_turnover_per_scenario <- function(w_matrix) {
  N <- nrow(w_matrix)
  T <- ncol(w_matrix)
  if (T < 2) return(rep(0, N))
  abs_diffs <- abs(w_matrix[, 2:T] - w_matrix[, 1:(T - 1)])
  rowSums(abs_diffs)
}

cat("Computing turnover for HybridNaive...\n")
turn_hybrid <- compute_turnover_per_scenario(w_hybrid)
cat(sprintf("  Turnover mean:   %.3f\n", mean(turn_hybrid)))
cat(sprintf("  Turnover median: %.3f\n\n", median(turn_hybrid)))

# ---- 6. Compute full metrics row for HybridNaive ------------------------
sharpe_like <- function(RR) {
  mean(RR) / sd(RR)
}

metrics_hybrid <- data.frame(
  strategy         = "HybridNaive",
  RR_mean          = round(mean(RR_hybrid), 3),
  RR_median        = round(median(RR_hybrid), 3),
  RR_sd            = round(sd(RR_hybrid), 3),
  RR_p10           = round(quantile(RR_hybrid, 0.10), 3),
  RR_p25           = round(quantile(RR_hybrid, 0.25), 3),
  RR_p75           = round(quantile(RR_hybrid, 0.75), 3),
  RR_p90           = round(quantile(RR_hybrid, 0.90), 3),
  RR_iqr           = round(quantile(RR_hybrid, 0.75) - quantile(RR_hybrid, 0.25), 3),
  P_shortfall_50pct = round(mean(RR_hybrid < 0.50), 3),
  P_shortfall_60pct = round(mean(RR_hybrid < 0.60), 3),
  P_shortfall_70pct = round(mean(RR_hybrid < 0.70), 3),
  P_shortfall_80pct = round(mean(RR_hybrid < 0.80), 3),
  CVaR_5pct        = round(mean(RR_hybrid[RR_hybrid <= quantile(RR_hybrid, 0.05)]), 3),
  CVaR_10pct       = round(mean(RR_hybrid[RR_hybrid <= quantile(RR_hybrid, 0.10)]), 3),
  RR_sharpe_like   = round(sharpe_like(RR_hybrid), 3),
  MDD_mean         = round(mean(mdd_hybrid), 3),
  MDD_median       = round(median(mdd_hybrid), 3),
  MDD_p10          = round(quantile(mdd_hybrid, 0.10), 3),
  turnover_mean    = round(mean(turn_hybrid), 3),
  turnover_median  = round(median(turn_hybrid), 3),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Strip the percentile names from quantile outputs
rownames(metrics_hybrid) <- NULL

# ---- 7. Append to existing structures -----------------------------------
metrics_new <- bind_rows(metrics_old, metrics_hybrid)

wp_new <- wp_old
wp_new[["HybridNaive"]] <- wealth_path_hybrid

mdd_new <- cbind(mdd_old, HybridNaive = mdd_hybrid)
to_new  <- cbind(to_old,  HybridNaive = turn_hybrid)
RR_new  <- cbind(RR_old,  HybridNaive = RR_hybrid)

# ---- 8. Display final comparison ----------------------------------------
cat("\n--- FULL METRICS TABLE (9 STRATEGIES) ---\n")
print(metrics_new %>%
        select(strategy, RR_median, RR_sd, P_shortfall_70pct,
               CVaR_5pct, MDD_mean, RR_sharpe_like, turnover_mean),
      row.names = FALSE)

# ---- 9. Save outputs (overwrites) ---------------------------------------
saveRDS(metrics_new, "Data/metrics_full.rds")
saveRDS(wp_new,      "Data/wealth_paths.rds")
saveRDS(mdd_new,     "Data/max_drawdowns.rds")
saveRDS(to_new,      "Data/turnovers.rds")
saveRDS(RR_new,      "Data/replacement_ratios_with_hybrid_naive.rds")

cat("\nSaved (overwrites with 9th strategy):\n")
cat("  Data/metrics_full.rds\n")
cat("  Data/wealth_paths.rds\n")
cat("  Data/max_drawdowns.rds\n")
cat("  Data/turnovers.rds\n")
cat("  Data/replacement_ratios_with_hybrid_naive.rds\n")
