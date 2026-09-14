library(dplyr)

cat("Loading inputs...\n")
metrics_old <- readRDS("Data/metrics_full.rds")
wp_old      <- readRDS("Data/wealth_paths.rds")
mdd_old     <- readRDS("Data/max_drawdowns.rds")
to_old      <- readRDS("Data/turnovers.rds")
RR_old      <- readRDS("Data/replacement_ratios_with_hybrid_naive.rds")

w_direct    <- readRDS("Data/direct_ml_allocations.rds")  # 20000 x 42, equity weights
eq_ret      <- readRDS("Data/equity_returns.rds")
bd_ret      <- readRDS("Data/bond_returns_10y.rds")
salary      <- readRDS("Data/salary_path.rds")
af          <- readRDS("Data/annuity_factor.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

cat("Direct_PPP equity weight stats (from saved allocations):\n")
cat(sprintf("  min:    %.3f\n", min(w_direct)))
cat(sprintf("  max:    %.3f\n", max(w_direct)))
cat(sprintf("  mean:   %.3f\n", mean(w_direct)))
cat(sprintf("  median: %.3f\n\n", median(w_direct)))

# Idempotency: strip previous Direct_PPP rows if present
if ("Direct_PPP" %in% metrics_old$strategy) {
  cat("Removing previous Direct_PPP row\n")
  metrics_old <- metrics_old %>% filter(strategy != "Direct_PPP")
  wp_old[["Direct_PPP"]] <- NULL
  mdd_old <- mdd_old[, colnames(mdd_old) != "Direct_PPP", drop = FALSE]
  to_old  <- to_old [, colnames(to_old)  != "Direct_PPP", drop = FALSE]
  RR_old  <- RR_old [, colnames(RR_old)  != "Direct_PPP", drop = FALSE]
}

# Wealth path
cat("Computing wealth path...\n")
wp_direct <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)
for (t in 1:N_YEARS) {
  contrib  <- CONTRIB_RATE * salary[, t]
  port_ret <- w_direct[, t] * eq_ret[, t] + (1 - w_direct[, t]) * bd_ret[, t]
  wp_direct[, t + 1] <- (wp_direct[, t] + contrib) * (1 + port_ret)
}
wf_direct <- wp_direct[, N_YEARS + 1]
RR_direct <- wf_direct / af / salary_at_retirement

# MDD
compute_mdd <- function(wp) {
  N <- nrow(wp); out <- numeric(N)
  for (i in 1:N) {
    path <- wp[i, ]; rm <- cummax(path)
    dd <- ifelse(rm > 0, (path - rm) / rm, 0)
    out[i] <- min(dd)
  }
  out
}
mdd_direct <- compute_mdd(wp_direct)

# Turnover
to_direct <- rowSums(abs(w_direct[, 2:N_YEARS] - w_direct[, 1:(N_YEARS - 1)]))

cat(sprintf("\nDirect_PPP stats (with NEW af):\n"))
cat(sprintf("  RR_median  = %.3f\n", median(RR_direct)))
cat(sprintf("  RR_mean    = %.3f\n", mean(RR_direct)))
cat(sprintf("  RR_sd      = %.3f\n", sd(RR_direct)))
cat(sprintf("  shortfall_70 = %.3f\n", mean(RR_direct < 0.70)))
cat(sprintf("  CVaR_5     = %.3f\n", mean(RR_direct[RR_direct <= quantile(RR_direct, 0.05)])))
cat(sprintf("  MDD_mean   = %.3f\n", mean(mdd_direct)))
cat(sprintf("  Turnover   = %.3f\n\n", mean(to_direct)))

# Build metrics row
metrics_direct <- data.frame(
  strategy = "Direct_PPP",
  RR_mean = round(mean(RR_direct), 3),
  RR_median = round(median(RR_direct), 3),
  RR_sd = round(sd(RR_direct), 3),
  RR_p10 = round(quantile(RR_direct, 0.10), 3),
  RR_p25 = round(quantile(RR_direct, 0.25), 3),
  RR_p75 = round(quantile(RR_direct, 0.75), 3),
  RR_p90 = round(quantile(RR_direct, 0.90), 3),
  RR_iqr = round(quantile(RR_direct, 0.75) - quantile(RR_direct, 0.25), 3),
  P_shortfall_50pct = round(mean(RR_direct < 0.50), 3),
  P_shortfall_60pct = round(mean(RR_direct < 0.60), 3),
  P_shortfall_70pct = round(mean(RR_direct < 0.70), 3),
  P_shortfall_80pct = round(mean(RR_direct < 0.80), 3),
  CVaR_5pct = round(mean(RR_direct[RR_direct <= quantile(RR_direct, 0.05)]), 3),
  CVaR_10pct = round(mean(RR_direct[RR_direct <= quantile(RR_direct, 0.10)]), 3),
  RR_sharpe_like = round(mean(RR_direct) / sd(RR_direct), 3),
  MDD_mean = round(mean(mdd_direct), 3),
  MDD_median = round(median(mdd_direct), 3),
  MDD_p10 = round(quantile(mdd_direct, 0.10), 3),
  turnover_mean = round(mean(to_direct), 3),
  turnover_median = round(median(to_direct), 3),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Append
metrics_new <- bind_rows(metrics_old, metrics_direct)
wp_new <- wp_old
wp_new[["Direct_PPP"]] <- wp_direct
mdd_new <- cbind(mdd_old, Direct_PPP = mdd_direct)
to_new  <- cbind(to_old, Direct_PPP = to_direct)
RR_new  <- cbind(RR_old, Direct_PPP = RR_direct)

cat("--- METRICS TABLE (10 STRATEGIES) ---\n")
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
cat("  Data/metrics_full.rds (10 strategies)\n")
