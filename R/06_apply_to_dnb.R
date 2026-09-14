# 06_apply_to_dnb.R  (V2 — fixes scale mismatch on yields and inflation)
#
# Input:  Data/ml_models.rds          (production Lasso, Ridge, RF)
#         Data/feature_cols.rds       (23 feature names in training order)
#         Data/mv_weights.rds         (static MV plain and LW weights)
#         Data/allocation_map.rds     (regime -> equity weight)
#         Data/equity_returns.rds     (20000 x 100)
#         Data/yield_1y.rds           (20000 x 101)
#         Data/yield_5y.rds           (20000 x 101)
#         Data/yield_10y.rds          (20000 x 101)
#         Data/yield_20y.rds          (20000 x 101)
#         Data/cp2022-...xlsx         (for EU inflation)
#
# Output: Data/dnb_equity_allocations.rds   (20000 x N_YEARS x 5 strategies)
#
# Dependencies: dplyr, readxl, glmnet, ranger
#
# ============================================================================
# IMPORTANT: SCALE CONVENTIONS
# ============================================================================
# Training (historical features in 02_feature_engineering.R) uses:
#   - Yields and inflation in PERCENTAGE POINTS (e.g. y10 = 2.35 meaning 2.35%)
#   - Equity returns in DECIMAL (e.g. eq_ret = 0.063 meaning +6.3%)
#   - eq_drawdown_12m in DECIMAL (e.g. -0.20 meaning -20%)
#   - dy1, dy10 in BASIS POINTS (Δyield × 100)
#
# DNB raw data (from Excel) uses DECIMAL throughout:
#   - hicp_yoy_dnb = 0.022 (meaning 2.2%)
#   - y10 = 0.0235 (meaning 2.35%)
#   - eq_ret = 0.063
#
# This script applies the conversion factor *100 to inflation and yields
# BEFORE building features, so that the input distribution matches training.
# Equity returns and drawdown are left in decimal (correct scale in both).
#
# Rule of thumb for downstream scripts: any feature that came out of HICP or
# yields needs to be in "pp" (percentage points). Equity-based features stay
# in decimal.
# ============================================================================

library(dplyr)
library(readxl)
library(glmnet)
library(ranger)

# ---- 1. Load all inputs ---------------------------------------------------
cat("Loading inputs...\n")

models       <- readRDS("Data/ml_models.rds")
feature_cols <- readRDS("Data/feature_cols.rds")
mv_weights   <- readRDS("Data/mv_weights.rds")
alloc_map    <- readRDS("Data/allocation_map.rds")

eq_ret  <- as.matrix(readRDS("Data/equity_returns.rds"))    # decimal
y1      <- as.matrix(readRDS("Data/yield_1y.rds"))          # decimal  -> rescale below
y5      <- as.matrix(readRDS("Data/yield_5y.rds"))          # decimal  -> rescale below
y10     <- as.matrix(readRDS("Data/yield_10y.rds"))         # decimal  -> rescale below
y20     <- as.matrix(readRDS("Data/yield_20y.rds"))         # decimal  -> rescale below

file_path <- "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx"
hicp_yoy_dnb <- as.matrix(read_excel(file_path,
                                     sheet = "5_EU_Price_Inflation",
                                     col_names = FALSE))     # decimal -> rescale below

# ---- 2. UNIT CONVERSION TO TRAINING SCALE -------------------------------
# Multiply yields and inflation by 100 to bring them from decimal to
# percentage points (training convention).
cat("Converting DNB yields and inflation from decimal to percentage points...\n")
y1           <- y1           * 100
y5           <- y5           * 100
y10          <- y10          * 100
y20          <- y20          * 100
hicp_yoy_dnb <- hicp_yoy_dnb * 100

cat("Post-conversion DNB scale check:\n")
cat(sprintf("  y10:          mean=%.3f  range=[%.3f, %.3f]\n",
            mean(y10), min(y10), max(y10)))
cat(sprintf("  hicp_yoy_dnb: mean=%.3f  range=[%.3f, %.3f]\n",
            mean(hicp_yoy_dnb), min(hicp_yoy_dnb), max(hicp_yoy_dnb)))
cat("Historical training scale (for reference):\n")
cat("  y10:          mean ~ 2.28  range=[-0.71,  5.35]\n")
cat("  hicp_yoy:     mean ~ 2.12  range=[-0.62, 10.09]\n\n")

# ---- 3. Configuration ----------------------------------------------------
N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42

cat("Strategy inference on", N_SCEN, "scenarios for", N_YEARS, "years\n\n")

# ---- 4. Feature builder --------------------------------------------------
# All inputs are now in TRAINING scale. No further conversion needed below.

expected_features <- c(
  "eq_ret", "hicp_yoy", "y1", "y5", "y10", "y20", "dy1", "dy10",
  "slope_10_1", "slope_20_5", "curvature", "real_yield10",
  "eq_drawdown_12m",
  "eq_ret_l1", "eq_ret_l3", "eq_ret_l12",
  "hicp_yoy_l3", "hicp_yoy_l12",
  "slope_l3", "slope_l12", "y10_l12",
  "slope_x_infl", "realy_x_eq"
)
stopifnot(identical(feature_cols, expected_features))

build_features_year <- function(t_career) {
  t_dnb_yield_curr <- t_career + 1
  t_dnb_yield_prev <- t_career
  t_dnb_ret_curr   <- t_career
  t_dnb_ret_prev   <- max(t_career - 1, 1)
  t_dnb_ret_l12    <- max(t_career - 1, 1)
  t_dnb_infl_curr  <- t_career
  t_dnb_infl_prev  <- max(t_career - 1, 1)
  
  # Current-year features (yields and inflation already in pp)
  eq_ret_t   <- eq_ret[, t_dnb_ret_curr]                  # decimal
  hicp_yoy_t <- hicp_yoy_dnb[, t_dnb_infl_curr]           # pp
  y1_t       <- y1[,  t_dnb_yield_curr]                   # pp
  y5_t       <- y5[,  t_dnb_yield_curr]
  y10_t      <- y10[, t_dnb_yield_curr]
  y20_t      <- y20[, t_dnb_yield_curr]
  
  # Yield changes in basis points
  if (t_career == 1) {
    y1_prev  <- y1[,  1]
    y10_prev <- y10[, 1]
  } else {
    y1_prev  <- y1[,  t_dnb_yield_prev]
    y10_prev <- y10[, t_dnb_yield_prev]
  }
  dy1_t  <- 100 * (y1_t  - y1_prev)
  dy10_t <- 100 * (y10_t - y10_prev)
  
  # Yield curve derived (all in pp)
  slope_10_1_t   <- y10_t - y1_t
  slope_20_5_t   <- y20_t - y5_t
  curvature_t    <- 2 * y5_t - y1_t - y10_t
  real_yield10_t <- y10_t - hicp_yoy_t
  
  # Equity drawdown — DECIMAL scale (matches training)
  eq_drawdown_12m_t <- pmin(0, eq_ret_t)
  
  # Lag features
  eq_ret_l1_t   <- if (t_career == 1) rep(0, N_SCEN) else eq_ret[, t_dnb_ret_prev]
  eq_ret_l3_t   <- eq_ret_l1_t
  eq_ret_l12_t  <- if (t_career == 1) rep(0, N_SCEN) else eq_ret[, t_dnb_ret_l12]
  hicp_yoy_l3_t  <- hicp_yoy_t
  hicp_yoy_l12_t <- if (t_career == 1) hicp_yoy_t else hicp_yoy_dnb[, t_dnb_infl_prev]
  
  if (t_career == 1) {
    slope_l3_t  <- slope_10_1_t
    slope_l12_t <- slope_10_1_t
    y10_l12_t   <- y10_t
  } else {
    y1_prev_full  <- y1[,  t_dnb_yield_prev]
    y10_prev_full <- y10[, t_dnb_yield_prev]
    slope_prev    <- y10_prev_full - y1_prev_full
    slope_l3_t    <- slope_prev
    slope_l12_t   <- slope_prev
    y10_l12_t     <- y10_prev_full
  }
  
  slope_x_infl_t <- slope_10_1_t * hicp_yoy_t
  realy_x_eq_t   <- real_yield10_t * eq_ret_l1_t
  
  X <- cbind(
    eq_ret           = eq_ret_t,
    hicp_yoy         = hicp_yoy_t,
    y1               = y1_t,
    y5               = y5_t,
    y10              = y10_t,
    y20              = y20_t,
    dy1              = dy1_t,
    dy10             = dy10_t,
    slope_10_1       = slope_10_1_t,
    slope_20_5       = slope_20_5_t,
    curvature        = curvature_t,
    real_yield10     = real_yield10_t,
    eq_drawdown_12m  = eq_drawdown_12m_t,
    eq_ret_l1        = eq_ret_l1_t,
    eq_ret_l3        = eq_ret_l3_t,
    eq_ret_l12       = eq_ret_l12_t,
    hicp_yoy_l3      = hicp_yoy_l3_t,
    hicp_yoy_l12     = hicp_yoy_l12_t,
    slope_l3         = slope_l3_t,
    slope_l12        = slope_l12_t,
    y10_l12          = y10_l12_t,
    slope_x_infl     = slope_x_infl_t,
    realy_x_eq       = realy_x_eq_t
  )
  stopifnot(identical(colnames(X), feature_cols))
  X
}

# ---- 5. Sanity check on year-5 features after rescaling -----------------
cat("--- POST-RESCALE CHECK: DNB FEATURES AT YEAR 5 ---\n")
X5 <- build_features_year(5)
cat(sprintf("hicp_yoy:    mean=%.3f  range=[%.3f, %.3f]\n",
            mean(X5[, "hicp_yoy"]), min(X5[, "hicp_yoy"]), max(X5[, "hicp_yoy"])))
cat(sprintf("slope_10_1:  mean=%.3f  range=[%.3f, %.3f]\n",
            mean(X5[, "slope_10_1"]), min(X5[, "slope_10_1"]), max(X5[, "slope_10_1"])))
cat(sprintf("y10:         mean=%.3f  range=[%.3f, %.3f]\n",
            mean(X5[, "y10"]), min(X5[, "y10"]), max(X5[, "y10"])))
cat("These should now be on the same scale as training (single-digit values).\n\n")

# ---- 6. Apply each strategy year by year --------------------------------
strategies <- c("MV_plain", "MV_LW", "Lasso", "Ridge", "RF")
allocations <- array(NA_real_,
                     dim = c(N_SCEN, N_YEARS, length(strategies)),
                     dimnames = list(NULL, NULL, strategies))

w_mv_plain <- as.numeric(mv_weights$plain["equity"])
w_mv_lw    <- as.numeric(mv_weights$lw["equity"])
allocations[, , "MV_plain"] <- w_mv_plain
allocations[, , "MV_LW"]    <- w_mv_lw

cat("MV_plain equity weight (static):", round(w_mv_plain, 3), "\n")
cat("MV_LW    equity weight (static):", round(w_mv_lw,    3), "\n\n")

cat("Predicting ML regimes year by year...\n")
start_time <- Sys.time()

for (t in 1:N_YEARS) {
  X_t <- build_features_year(t)
  
  prob_lasso  <- predict(models$lasso, newx = X_t,
                         s = "lambda.min", type = "response")[, , 1]
  class_lasso <- colnames(prob_lasso)[apply(prob_lasso, 1, which.max)]
  allocations[, t, "Lasso"] <- alloc_map[class_lasso]
  
  prob_ridge  <- predict(models$ridge, newx = X_t,
                         s = "lambda.min", type = "response")[, , 1]
  class_ridge <- colnames(prob_ridge)[apply(prob_ridge, 1, which.max)]
  allocations[, t, "Ridge"] <- alloc_map[class_ridge]
  
  prob_rf  <- predict(models$rf, data = as.data.frame(X_t))$predictions
  class_rf <- colnames(prob_rf)[apply(prob_rf, 1, which.max)]
  allocations[, t, "RF"] <- alloc_map[class_rf]
  
  if (t %% 5 == 0) {
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    cat(sprintf("  Year %2d / %d done (%.1fs)\n", t, N_YEARS, elapsed))
  }
}

elapsed_total <- as.numeric(Sys.time() - start_time, units = "secs")
cat(sprintf("\nML inference completed in %.1fs.\n\n", elapsed_total))

# ---- 7. Diagnostics ------------------------------------------------------
cat("--- ALLOCATION DIAGNOSTICS ---\n\n")

for (s in strategies) {
  cat(s, ":\n", sep = "")
  alloc_s <- allocations[, , s]
  cat(sprintf("  Mean:    %.3f\n", mean(alloc_s)))
  cat(sprintf("  Median:  %.3f\n", median(alloc_s)))
  cat(sprintf("  SD:      %.3f\n", sd(alloc_s)))
  cat(sprintf("  Range:   [%.3f, %.3f]\n", min(alloc_s), max(alloc_s)))
  cat("\n")
}

cat("--- ML EQUITY ALLOCATION DISTRIBUTION (% across scen-years) ---\n")
unique_alloc <- sort(unique(c(allocations[, , "Lasso"],
                              allocations[, , "Ridge"],
                              allocations[, , "RF"])))
dist_table <- sapply(c("Lasso", "Ridge", "RF"), function(s) {
  tab <- table(factor(round(allocations[, , s], 2),
                      levels = round(unique_alloc, 2)))
  round(prop.table(tab) * 100, 1)
})
print(dist_table)

cat("\n--- MEAN EQUITY ALLOCATION BY YEAR (first 10 years) ---\n")
yearly_means <- sapply(c("Lasso", "Ridge", "RF"), function(s) {
  round(colMeans(allocations[, 1:10, s]), 3)
})
print(yearly_means)

# ---- 8. Save -------------------------------------------------------------
saveRDS(allocations, "Data/dnb_equity_allocations.rds")
cat("\nSaved Data/dnb_equity_allocations.rds\n")
cat("Dimension:", dim(allocations), "\n")

cat("--- AG2024_cohort_survival ---\n")
surv <- readRDS("Data/AG2024_cohort_survival.rds")
cat("Class:", class(surv), "\n")
cat("Length/dim:", if(is.null(dim(surv))) length(surv) else dim(surv), "\n")
cat("Range:", round(range(surv), 4), "\n")

cat("\n--- annuity_factor ---\n")
af <- readRDS("Data/annuity_factor.rds")
cat("Class:", class(af), "\n")
cat("Dim:", if(is.null(dim(af))) length(af) else dim(af), "\n")
cat("Range:", round(range(af, na.rm=TRUE), 4), "\n")
if(!is.null(dim(af))) {
  cat("First 5 values:", round(head(as.vector(af), 5), 4), "\n")
} else {
  cat("Value(s):", round(head(af, 5), 4), "\n")
}

cat("\n--- salary_path ---\n")
sp <- readRDS("Data/salary_path.rds")
cat("Class:", class(sp), "\n")
cat("Dim:", dim(sp), "\n")
cat("Range:", round(range(sp), 4), "\n")
cat("First 3 scenarios, first 5 years:\n")
print(round(sp[1:3, 1:5], 2))
cat("First 3 scenarios, years 40-42:\n")
print(round(sp[1:3, 40:42], 2))

cat("\n--- replacement_ratios_modelfree ---\n")
rr <- readRDS("Data/replacement_ratios_modelfree.rds")
cat("Class:", class(rr), "\n")
cat("Dim:", if(is.null(dim(rr))) length(rr) else dim(rr), "\n")
if (is.list(rr)) {
  cat("Names:", names(rr), "\n")
  for (n in names(rr)) {
    cat("  ", n, ": class=", class(rr[[n]]), " len/dim=",
        if(is.null(dim(rr[[n]]))) length(rr[[n]]) else paste(dim(rr[[n]]), collapse="x"),
        "\n", sep="")
  }
} else if (is.matrix(rr)) {
  cat("First 3 rows, first 5 cols:\n")
  print(round(rr[1:3, 1:5], 4))
  cat("Median by column (strategy):\n")
  print(round(apply(rr, 2, median), 4))
}

cat("--- AG2024_cohort_survival ---\n")
surv <- readRDS("Data/AG2024_cohort_survival.rds")
cat("Class:", class(surv), "\n")
cat("Structure:\n"); str(surv)
cat("\nHead (raw):\n"); print(head(surv, 10))
cat("\nTail (raw):\n"); print(tail(surv, 10))

cat("\n--- annuity_factor ---\n")
af <- readRDS("Data/annuity_factor.rds")
cat("Class:", class(af), "\n")
cat("Structure:\n"); str(af)
if (is.numeric(af) && is.null(dim(af))) {
  cat("Range:", round(range(af, na.rm=TRUE), 4), "\n")
  cat("Head:", round(head(af, 5), 4), "\n")
}

cat("\n--- salary_path ---\n")
sp <- readRDS("Data/salary_path.rds")
cat("Class:", class(sp), "\n")
cat("Dim:", dim(sp), "\n")
cat("Structure:\n"); str(sp, max.level = 1)
if (is.numeric(sp)) {
  cat("Range:", round(range(sp, na.rm=TRUE), 2), "\n")
  cat("First 3 rows, first 5 cols:\n")
  print(round(sp[1:3, 1:5], 2))
  cat("First 3 rows, last 3 cols:\n")
  print(round(sp[1:3, (ncol(sp)-2):ncol(sp)], 2))
}

cat("\n--- replacement_ratios_modelfree ---\n")
rr <- readRDS("Data/replacement_ratios_modelfree.rds")
cat("Class:", class(rr), "\n")
cat("Structure:\n"); str(rr, max.level = 1)
