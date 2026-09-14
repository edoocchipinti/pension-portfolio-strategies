library(dplyr)

cat("Loading inputs...\n")
metrics_old <- readRDS("Data/metrics_full.rds")
wp_old      <- readRDS("Data/wealth_paths.rds")
mdd_old     <- readRDS("Data/max_drawdowns.rds")
to_old      <- readRDS("Data/turnovers.rds")
RR_old      <- readRDS("Data/replacement_ratios_all.rds")

allocations <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret      <- readRDS("Data/equity_returns.rds")
bd_ret      <- readRDS("Data/bond_returns_10y.rds")
salary      <- readRDS("Data/salary_path.rds")
af          <- readRDS("Data/annuity_factor.rds")
mv_weights  <- readRDS("Data/mv_weights.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

cat("Existing strategies in metrics_full:\n")
print(metrics_old$strategy)
cat("Adding: HybridNaive\n\n")

# Build HybridNaive equity weight matrix: (1/3) * 60/40 + (1/3) * MV_plain + (1/3) * RF
w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_weights$plain["equity"]), N_SCEN, N_YEARS)
w_rf    <- allocations[, , "RF"]
w_hybrid <- (1/3) * w_6040 + (1/3) * w_mv_p + (1/3) * w_rf

cat("HybridNaive equity weight stats:\n")
cat(sprintf("  min:    %.3f\n", min(w_hybrid)))
cat(sprintf("  max:    %.3f\n", max(w_hybrid)))
cat(sprintf("  mean:   %.3f\n\n", mean(w_hybrid)))

# Wealth path
wp_hybrid <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)
for (t in 1:N_YEARS) {
  contrib  <- CONTRIB_RATE * salary[, t]
  port_ret <- w_hybrid[, t] * eq_ret[, t] + (1 - w_hybrid[, t]) * bd_ret[, t]
  wp_hybrid[, t + 1] <- (wp_hybrid[, t] + contrib) * (1 + port_ret)
}
wf_hybrid <- wp_hybrid[, N_YEARS + 1]
RR_hybrid <- wf_hybrid / af / salary_at_retirement

cat("HybridNaive RR stats:\n")
cat(sprintf("  mean:   %.3f\n", mean(RR_hybrid)))
cat(sprintf("  median: %.3f\n", median(RR_hybrid)))
cat(sprintf("  sd:     %.3f\n\n", sd(RR_hybrid)))

# MDD
compute_mdd <- function(wp) {
  N <- nrow(wp)
  mdd <- numeric(N)
  for (i in 1:N) {
    path <- wp[i, ]
    rm <- cummax(path)
    dd <- ifelse(rm > 0, (path - rm) / rm, 0)
    mdd[i] <- min(dd)
  }
  mdd
}
mdd_hybrid <- compute_mdd(wp_hybrid)

# Turnover
delta <- abs(w_hybrid[, 2:N_YEARS] - w_hybrid[, 1:(N_YEARS - 1)])
to_hybrid <- rowSums(delta)

cat(sprintf("MDD mean:      %.3f\n", mean(mdd_hybrid)))
cat(sprintf("Turnover mean: %.3f\n\n", mean(to_hybrid)))

# Build metrics row
metrics_hybrid <- data.frame(
  strategy = "HybridNaive",
  RR_mean = round(mean(RR_hybrid), 3),
  RR_median = round(median(RR_hybrid), 3),
  RR_sd = round(sd(RR_hybrid), 3),
  RR_p10 = round(quantile(RR_hybrid, 0.10), 3),
  RR_p25 = round(quantile(RR_hybrid, 0.25), 3),
  RR_p75 = round(quantile(RR_hybrid, 0.75), 3),
  RR_p90 = round(quantile(RR_hybrid, 0.90), 3),
  RR_iqr = round(quantile(RR_hybrid, 0.75) - quantile(RR_hybrid, 0.25), 3),
  P_shortfall_50pct = round(mean(RR_hybrid < 0.50), 3),
  P_shortfall_60pct = round(mean(RR_hybrid < 0.60), 3),
  P_shortfall_70pct = round(mean(RR_hybrid < 0.70), 3),
  P_shortfall_80pct = round(mean(RR_hybrid < 0.80), 3),
  CVaR_5pct = round(mean(RR_hybrid[RR_hybrid <= quantile(RR_hybrid, 0.05)]), 3),
  CVaR_10pct = round(mean(RR_hybrid[RR_hybrid <= quantile(RR_hybrid, 0.10)]), 3),
  RR_sharpe_like = round(mean(RR_hybrid) / sd(RR_hybrid), 3),
  MDD_mean = round(mean(mdd_hybrid), 3),
  MDD_median = round(median(mdd_hybrid), 3),
  MDD_p10 = round(quantile(mdd_hybrid, 0.10), 3),
  turnover_mean = round(mean(to_hybrid), 3),
  turnover_median = round(median(to_hybrid), 3),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Append
metrics_new <- bind_rows(metrics_old, metrics_hybrid)
wp_new <- wp_old
wp_new[["HybridNaive"]] <- wp_hybrid
mdd_new <- cbind(mdd_old, HybridNaive = mdd_hybrid)
to_new  <- cbind(to_old, HybridNaive = to_hybrid)
RR_new  <- cbind(RR_old, HybridNaive = RR_hybrid)

cat("--- METRICS TABLE (9 STRATEGIES) ---\n")
print(metrics_new[, c("strategy", "RR_median", "RR_sd", "P_shortfall_70pct", 
                      "CVaR_5pct", "MDD_mean", "turnover_mean")],
      row.names = FALSE)

# Save
saveRDS(metrics_new, "Data/metrics_full.rds")
saveRDS(wp_new, "Data/wealth_paths.rds")
saveRDS(mdd_new, "Data/max_drawdowns.rds")
saveRDS(to_new, "Data/turnovers.rds")
saveRDS(RR_new, "Data/replacement_ratios_with_hybrid_naive.rds")

cat("\n--- SAVED ---\n")
cat("  Data/metrics_full.rds (9 strategies)\n")
cat("  Data/wealth_paths.rds\n")
cat("  Data/max_drawdowns.rds\n")
cat("  Data/turnovers.rds\n")
cat("  Data/replacement_ratios_with_hybrid_naive.rds\n")
