# 07_compute_rr.R
#
# Input:  Data/dnb_equity_allocations.rds       (20000 x 42 x 5)
#         Data/equity_returns.rds               (20000 x 100)
#         Data/bond_returns_10y.rds             (20000 x 100)
#         Data/salary_path.rds                  (20000 x 43; t=0 to t=42)
#         Data/annuity_factor.rds               (20000)
#         Data/AG2024_cohort_survival.rds       (list)
#         Data/replacement_ratios_modelfree.rds (list with RR_6040, RR_1N, RR_glide)
#
# Output: Data/replacement_ratios_all.rds       (matrix 20000 x 8 strategies)
#
# Dependencies: none beyond base R
#
# Computes the replacement ratio (RR) for all 8 strategies on the 20,000 DNB
# scenarios. The RR is the ratio between annual pension annuity and final
# salary at retirement.
#
# Strategies:
#   Rule-based (static / age-dependent):  60/40, 1/N, glide-path
#   Statistical optimisation (static):    MV plain, MV LW
#   Machine learning (regime-dynamic):    Lasso, Ridge, RF
#
# Cross-check: 60/40 RR computed here should match RR_6040 from Phase 2.

# ---- 1. Load inputs -------------------------------------------------------
cat("Loading inputs...\n")
allocations <- readRDS("Data/dnb_equity_allocations.rds")  # 20000 x 42 x 5
eq_ret      <- readRDS("Data/equity_returns.rds")          # 20000 x 100
bd_ret      <- readRDS("Data/bond_returns_10y.rds")        # 20000 x 100
salary      <- readRDS("Data/salary_path.rds")             # 20000 x 43
af          <- readRDS("Data/annuity_factor.rds")          # 20000
rr_phase2   <- readRDS("Data/replacement_ratios_modelfree.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]   # 42

cat(sprintf("N_SCEN = %d, N_YEARS = %d\n", N_SCEN, N_YEARS))

# ---- 2. Define ALL 8 strategies' equity weight paths ---------------------
# Returns a list of 20000 x 42 matrices, one per strategy.

# Rule-based: build from scratch (same logic as Phase 2)
w_6040  <- matrix(0.60, nrow = N_SCEN, ncol = N_YEARS)
w_1N    <- matrix(0.50, nrow = N_SCEN, ncol = N_YEARS)
# Glide-path: linear from 0.80 at age 25 (year 1) to 0.30 at age 67 (year 42)
glide_path_vec <- seq(0.80, 0.30, length.out = N_YEARS)
w_glide <- matrix(rep(glide_path_vec, each = N_SCEN), nrow = N_SCEN, ncol = N_YEARS)

# MV (static)
w_mv_plain <- allocations[, , "MV_plain"]
w_mv_lw    <- allocations[, , "MV_LW"]

# ML (regime-dynamic)
w_lasso <- allocations[, , "Lasso"]
w_ridge <- allocations[, , "Ridge"]
w_rf    <- allocations[, , "RF"]

strategy_weights <- list(
  "60/40"     = w_6040,
  "1/N"       = w_1N,
  "Glide"     = w_glide,
  "MV_plain"  = w_mv_plain,
  "MV_LW"     = w_mv_lw,
  "Lasso"     = w_lasso,
  "Ridge"     = w_ridge,
  "RF"        = w_rf
)

cat("\nMean equity weight per strategy (across all scenario-years):\n")
for (s in names(strategy_weights)) {
  cat(sprintf("  %-10s  %.3f\n", s, mean(strategy_weights[[s]])))
}

# ---- 3. Wealth accumulation function -------------------------------------
# For each scenario i and strategy s:
#   contribution_t = 0.14 * salary[i, t]   (t = 1..42 mapping to indices 1..42)
#   portfolio_ret_t = w_eq[i,t] * eq_ret[i,t] + (1-w_eq[i,t]) * bond_ret[i,t]
#   wealth[i, t] = (wealth[i, t-1] + contribution_t) * (1 + portfolio_ret_t)
# Final: W_T = wealth[i, 42]
#
# Convention check: salary_path[, 1] = 1 at age 25.
# The "current" salary in year t (during year t of career) is salary[, t].
# Contribution is paid at the start of year t, then invested at portfolio_ret_t.

CONTRIB_RATE <- 0.14

compute_wealth_final <- function(w_eq) {
  # w_eq: N_SCEN x N_YEARS
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 4. Compute final wealth per strategy --------------------------------
cat("\nComputing terminal wealth for each strategy...\n")
final_wealth <- sapply(strategy_weights, compute_wealth_final)
# final_wealth: N_SCEN x 8

cat("Mean terminal wealth (relative to initial salary = 1):\n")
print(round(colMeans(final_wealth), 3))
cat("Median terminal wealth:\n")
print(round(apply(final_wealth, 2, median), 3))

# ---- 5. Compute Replacement Ratio ----------------------------------------
# RR[i, s] = (W_T[i, s] / annuity_factor[i]) / salary_at_retirement[i]
# salary_at_retirement = salary[, 43] (after 42 years of growth)

salary_at_retirement <- salary[, N_YEARS + 1]  # index 43 = age 67 salary

annual_annuity <- sweep(final_wealth, 1, af, FUN = "/")
RR <- sweep(annual_annuity, 1, salary_at_retirement, FUN = "/")
colnames(RR) <- names(strategy_weights)

# ---- 6. CROSS-CHECK against Phase 2 --------------------------------------
cat("\n--- CROSS-CHECK AGAINST PHASE 2 ---\n")
cat("This verifies that the new wealth-accumulation formula matches Phase 2.\n\n")

# Compare 60/40
cor_6040 <- cor(RR[, "60/40"], rr_phase2$RR_6040)
diff_6040 <- RR[, "60/40"] - rr_phase2$RR_6040
cat(sprintf("60/40:  median new=%.4f  Phase2=%.4f  correlation=%.4f  max_abs_diff=%.4f\n",
            median(RR[, "60/40"]), median(rr_phase2$RR_6040),
            cor_6040, max(abs(diff_6040))))

cor_1N <- cor(RR[, "1/N"], rr_phase2$RR_1N)
diff_1N <- RR[, "1/N"] - rr_phase2$RR_1N
cat(sprintf("1/N:    median new=%.4f  Phase2=%.4f  correlation=%.4f  max_abs_diff=%.4f\n",
            median(RR[, "1/N"]), median(rr_phase2$RR_1N),
            cor_1N, max(abs(diff_1N))))

cor_glide <- cor(RR[, "Glide"], rr_phase2$RR_glide)
diff_glide <- RR[, "Glide"] - rr_phase2$RR_glide
cat(sprintf("Glide:  median new=%.4f  Phase2=%.4f  correlation=%.4f  max_abs_diff=%.4f\n",
            median(RR[, "Glide"]), median(rr_phase2$RR_glide),
            cor_glide, max(abs(diff_glide))))

cat("\nExpected from handover: 60/40 median=0.5447, 1/N=0.5013, Glide=0.4857\n")

# ---- 7. Headline metrics for all 8 strategies ----------------------------
cat("\n--- REPLACEMENT RATIO BY STRATEGY ---\n")

RR_TARGET <- 0.70

summary_table <- data.frame(
  strategy        = colnames(RR),
  mean            = round(colMeans(RR),               3),
  median          = round(apply(RR, 2, median),        3),
  sd              = round(apply(RR, 2, sd),            3),
  p10             = round(apply(RR, 2, quantile, 0.10), 3),
  p25             = round(apply(RR, 2, quantile, 0.25), 3),
  p75             = round(apply(RR, 2, quantile, 0.75), 3),
  p90             = round(apply(RR, 2, quantile, 0.90), 3),
  shortfall_70pct = round(colMeans(RR < RR_TARGET),    3)
)
print(summary_table, row.names = FALSE)

# ---- 8. Save -------------------------------------------------------------
saveRDS(RR, "Data/replacement_ratios_all.rds")
saveRDS(summary_table, "Data/rr_summary_table.rds")
cat("\nSaved Data/replacement_ratios_all.rds\n")
cat("Saved Data/rr_summary_table.rds\n")
