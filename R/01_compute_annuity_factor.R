# 01_compute_annuity_factor.R
#
# Input:  Data/AG2024_cohort_survival.rds   (cohort 1999 survival probabilities)
#         Data/yield_1y.rds                  (20000 x 101, decimal)
#         Data/yield_5y.rds                  (20000 x 101, decimal)
#         Data/yield_10y.rds                 (20000 x 101, decimal)
#         Data/yield_20y.rds                 (20000 x 101, decimal)
#
# Output: Data/annuity_factor_v2.rds         (20000-vector, scenario-specific)
#         Data/annuity_factor_diagnostic.rds (comparison with af_v1)
#
# Dependencies: none beyond base R
#
# Computes the scenario-specific annuity factor for all 20,000 DNB
# scenarios at retirement (year 42, age 67) using the standard actuarial
# present-value formula:
#
#   a_i = sum_{k=0}^{K-1}  P(T_death > 67+k) / (1 + y_{k,i})^k
#
# where:
#   - P(T_death > 67+k): survival probability from AG2024 cohort 1999 table
#     (k=0 corresponds to age 67, k=K-1 to age 67+K-1 = 111)
#   - y_{k,i}: zero-coupon yield at maturity k years, scenario i, year 42
#     (linearly interpolated from yields at maturities 1, 5, 10, 20)
#   - K = 45 (annuity payable from age 67 to age 111)
#
# Conventions:
#   - Annuity-due (first payment at retirement, k=0)
#   - Yields in decimal scale (e.g., 0.0282 = 2.82%)
#   - Linear interpolation of yield curve between {1,5,10,20} year maturities
#     For k = 0: use y1 (short-end extrapolation)
#     For k > 20: use y20 (long-end extrapolation, matches DNB asymptotic
#                          term structure behavior)
#
# Output is saved as annuity_factor_v2.rds to preserve the original
# (annuity_factor.rds) for diagnostic comparison. To switch the pipeline
# to v2, manually rename or update 07_compute_rr.R to read v2.

# ---- 1. Configuration ---------------------------------------------------
N_SCEN <- 20000
RETIREMENT_YEAR_COL <- 42   # column index in yield matrices at retirement
K <- 45                      # number of annuity payment years (ages 67..111)
YIELD_MATURITIES <- c(1, 5, 10, 20)

cat("================================================================\n")
cat("ANNUITY FACTOR COMPUTATION (V2 - standard actuarial formula)\n")
cat("================================================================\n\n")

cat("Configuration:\n")
cat(sprintf("  N scenarios:           %d\n",        N_SCEN))
cat(sprintf("  Retirement year col:   %d\n",        RETIREMENT_YEAR_COL))
cat(sprintf("  Annuity payment years: K = %d (ages 67 to %d)\n",
            K, 67 + K - 1))
cat(sprintf("  Yield maturities:      {%s}\n\n",
            paste(YIELD_MATURITIES, collapse = ", ")))

# ---- 2. Load inputs ------------------------------------------------------
cat("Loading inputs...\n")
cohort  <- readRDS("Data/AG2024_cohort_survival.rds")
y1      <- as.matrix(readRDS("Data/yield_1y.rds"))
y5      <- as.matrix(readRDS("Data/yield_5y.rds"))
y10     <- as.matrix(readRDS("Data/yield_10y.rds"))
y20     <- as.matrix(readRDS("Data/yield_20y.rds"))
af_v1   <- readRDS("Data/annuity_factor.rds")  # for diagnostic comparison

survival_probs <- cohort$survival
stopifnot(length(survival_probs) >= K)

cat(sprintf("  Survival probabilities: length=%d (need >= %d)\n",
            length(survival_probs), K))
cat(sprintf("  Yield matrices:         %d x %d each\n",
            nrow(y10), ncol(y10)))
cat(sprintf("  Cohort retirement year: %d\n",   cohort$retirement_year))
cat(sprintf("  Cohort e67:             %.4f years\n\n", cohort$e67))

# ---- 3. Yield curve interpolation function ------------------------------
# Returns the interpolated yield at maturity tau (in years), given
# the yield curve values at maturities {1, 5, 10, 20}.
# - For tau < 1: returns y1 (short-end flat extrapolation)
# - For tau in [1, 20]: linear interpolation
# - For tau > 20: returns y20 (long-end flat extrapolation)
#
# Vectorised version: works on full vector of scenarios at once.
interpolate_yield_vec <- function(tau, y1_v, y5_v, y10_v, y20_v) {
  if (tau <= 1)  return(y1_v)
  if (tau >= 20) return(y20_v)
  
  if (tau <= 5) {
    # Linear between y1 and y5
    w <- (tau - 1) / (5 - 1)
    return((1 - w) * y1_v + w * y5_v)
  }
  if (tau <= 10) {
    # Linear between y5 and y10
    w <- (tau - 5) / (10 - 5)
    return((1 - w) * y5_v + w * y10_v)
  }
  # tau in (10, 20)
  w <- (tau - 10) / (20 - 10)
  (1 - w) * y10_v + w * y20_v
}

# ---- 4. Compute annuity factor for all scenarios ------------------------
# Vectorised computation: for each k from 0 to K-1, compute the discount
# factor across all scenarios at once, then accumulate the contribution
# survival[k+1] / (1 + y_k)^k

cat("Computing annuity factor for all", N_SCEN, "scenarios...\n")
start_time <- Sys.time()

# Extract yields at retirement (year 42) for all scenarios
y1_ret  <- y1[,  RETIREMENT_YEAR_COL]
y5_ret  <- y5[,  RETIREMENT_YEAR_COL]
y10_ret <- y10[, RETIREMENT_YEAR_COL]
y20_ret <- y20[, RETIREMENT_YEAR_COL]

# Initialise accumulator
af_v2 <- numeric(N_SCEN)

for (k in 0:(K - 1)) {
  # Yield at maturity k, for all scenarios at once
  y_k <- interpolate_yield_vec(k, y1_ret, y5_ret, y10_ret, y20_ret)
  
  # Survival probability at age 67+k (k=0 -> survival[1] = 1)
  surv <- survival_probs[k + 1]
  
  # Discount factor: (1 + y_k)^k (k=0 -> discount = 1, payable immediately)
  discount <- (1 + y_k)^k
  
  # Accumulate
  af_v2 <- af_v2 + surv / discount
}

elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
cat(sprintf("Computation completed in %.2fs\n\n", elapsed))

# ---- 5. Sanity checks ---------------------------------------------------
cat("--- ANNUITY FACTOR DISTRIBUTION (v2 - new) ---\n")
cat(sprintf("  Min:     %.4f\n", min(af_v2)))
cat(sprintf("  Q1:      %.4f\n", quantile(af_v2, 0.25)))
cat(sprintf("  Median:  %.4f\n", median(af_v2)))
cat(sprintf("  Mean:    %.4f\n", mean(af_v2)))
cat(sprintf("  Q3:      %.4f\n", quantile(af_v2, 0.75)))
cat(sprintf("  Max:     %.4f\n", max(af_v2)))
cat(sprintf("  SD:      %.4f\n\n", sd(af_v2)))

cat("--- ANNUITY FACTOR DISTRIBUTION (v1 - existing, for comparison) ---\n")
cat(sprintf("  Min:     %.4f\n", min(af_v1)))
cat(sprintf("  Q1:      %.4f\n", quantile(af_v1, 0.25)))
cat(sprintf("  Median:  %.4f\n", median(af_v1)))
cat(sprintf("  Mean:    %.4f\n", mean(af_v1)))
cat(sprintf("  Q3:      %.4f\n", quantile(af_v1, 0.75)))
cat(sprintf("  Max:     %.4f\n", max(af_v1)))
cat(sprintf("  SD:      %.4f\n\n", sd(af_v1)))

# ---- 6. Compare v1 vs v2 ------------------------------------------------
cat("--- V1 vs V2 COMPARISON ---\n")
diff_abs <- af_v2 - af_v1
diff_rel <- diff_abs / af_v1

cat(sprintf("  Mean abs diff:    %.4f\n",   mean(diff_abs)))
cat(sprintf("  Median abs diff:  %.4f\n",   median(diff_abs)))
cat(sprintf("  Mean rel diff:    %.2f%%\n", 100 * mean(diff_rel)))
cat(sprintf("  Median rel diff:  %.2f%%\n", 100 * median(diff_rel)))
cat(sprintf("  Correlation:      %.4f\n",   cor(af_v1, af_v2)))

cat("\nScatter check (first 5 scenarios):\n")
comparison_head <- data.frame(
  scenario = 1:5,
  af_v1    = round(af_v1[1:5], 4),
  af_v2    = round(af_v2[1:5], 4),
  diff     = round(diff_abs[1:5], 4),
  rel_diff = round(100 * diff_rel[1:5], 2)
)
print(comparison_head, row.names = FALSE)

# ---- 7. Interpretation guidance -----------------------------------------
mean_rel <- mean(diff_rel)
cat("\n--- INTERPRETATION ---\n")
if (abs(mean_rel) < 0.01) {
  cat("V2 essentially matches V1 (mean rel diff < 1%). Safe to switch.\n")
} else if (abs(mean_rel) < 0.05) {
  cat("V2 differs slightly from V1 (mean rel diff < 5%). Investigate but\n")
  cat("V2 is likely correct (standard textbook formula).\n")
} else {
  cat("V2 differs notably from V1 (mean rel diff >= 5%). V2 follows the\n")
  cat("standard actuarial formula (Bowers et al. 1997; Bühlmann 1992).\n")
  cat("Difference suggests V1 used a non-standard convention.\n")
}

if (cor(af_v1, af_v2) > 0.95) {
  cat("Correlation > 0.95: same underlying signal, only scaling differs.\n")
  cat("Rankings across scenarios will be preserved when switching to V2.\n")
} else if (cor(af_v1, af_v2) > 0.80) {
  cat("Correlation 0.80-0.95: similar signal but some scenarios will\n")
  cat("re-rank when switching to V2.\n")
} else {
  cat("Correlation < 0.80: V1 and V2 disagree structurally. Investigate.\n")
}

# ---- 8. Save outputs ----------------------------------------------------
saveRDS(af_v2, "Data/annuity_factor_v2.rds")
cat("\nSaved Data/annuity_factor_v2.rds (20,000-vector)\n")

# Diagnostic file for thesis Appendix C
diagnostic <- list(
  af_v2          = af_v2,
  af_v1          = af_v1,
  diff_abs       = diff_abs,
  diff_rel       = diff_rel,
  correlation    = cor(af_v1, af_v2),
  config         = list(
    N_SCEN              = N_SCEN,
    RETIREMENT_YEAR_COL = RETIREMENT_YEAR_COL,
    K                   = K,
    YIELD_MATURITIES    = YIELD_MATURITIES,
    formula             = "a = sum_{k=0}^{K-1} P(T>67+k) / (1+y_k)^k",
    convention          = "annuity-due, end-of-year yield, linear interpolation"
  )
)
saveRDS(diagnostic, "Data/annuity_factor_diagnostic.rds")
cat("Saved Data/annuity_factor_diagnostic.rds (for Appendix C)\n\n")

cat("================================================================\n")
cat("NEXT STEPS:\n")
cat("  1. Review the v1 vs v2 comparison above.\n")
cat("  2. If satisfied, modify R/07_compute_rr.R to load v2:\n")
cat("       af <- readRDS('Data/annuity_factor_v2.rds')\n")
cat("  3. Rerun pipeline:\n")
cat("       source('R/07_compute_rr.R')\n")
cat("       source('R/08_metrics.R')\n")
cat("       source('R/09_plots.R')\n")
cat("       source('R/11_hybrid_strategy.R')\n")
cat("       (and any robustness scripts that use af)\n")
cat("================================================================\n")
