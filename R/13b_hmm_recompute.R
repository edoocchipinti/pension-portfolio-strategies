library(dplyr)

cat("Loading inputs...\n")
metrics_old <- readRDS("Data/metrics_full.rds")
wp_old      <- readRDS("Data/wealth_paths.rds")
mdd_old     <- readRDS("Data/max_drawdowns.rds")
to_old      <- readRDS("Data/turnovers.rds")
RR_old      <- readRDS("Data/replacement_ratios_with_hybrid_naive.rds")

hmm_alloc   <- readRDS("Data/hmm_dnb_allocations.rds")  # list: hard, soft, states
eq_ret      <- readRDS("Data/equity_returns.rds")
bd_ret      <- readRDS("Data/bond_returns_10y.rds")
salary      <- readRDS("Data/salary_path.rds")
af          <- readRDS("Data/annuity_factor.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

w_hard <- hmm_alloc$hard
w_soft <- hmm_alloc$soft

cat("HMM_hard stats: mean=", round(mean(w_hard), 3), 
    " range=[", round(min(w_hard), 3), ",", round(max(w_hard), 3), "]\n")
cat("HMM_soft stats: mean=", round(mean(w_soft), 3),
    " range=[", round(min(w_soft), 3), ",", round(max(w_soft), 3), "]\n\n")

# Idempotency: strip duplicates AND previous HMM rows
metrics_old <- metrics_old %>% 
  filter(!duplicated(strategy)) %>%
  filter(!strategy %in% c("HMM_hard", "HMM_soft"))

# Also remove from other structures
for (s in c("HMM_hard", "HMM_soft")) {
  wp_old[[s]] <- NULL
  mdd_old <- mdd_old[, colnames(mdd_old) != s, drop = FALSE]
  to_old  <- to_old [, colnames(to_old)  != s, drop = FALSE]
  RR_old  <- RR_old [, colnames(RR_old)  != s, drop = FALSE]
}

cat("Current metrics_full strategies after dedup:\n")
print(metrics_old$strategy)
cat("\n")

# Compute wealth, MDD, turnover for both HMM
compute_wealth_path <- function(w_eq) {
  wp <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wp[, t + 1] <- (wp[, t] + contrib) * (1 + port_ret)
  }
  wp
}

compute_mdd <- function(wp) {
  N <- nrow(wp); out <- numeric(N)
  for (i in 1:N) {
    path <- wp[i, ]; rm <- cummax(path)
    dd <- ifelse(rm > 0, (path - rm) / rm, 0)
    out[i] <- min(dd)
  }
  out
}

compute_turnover <- function(w_eq) {
  rowSums(abs(w_eq[, 2:ncol(w_eq)] - w_eq[, 1:(ncol(w_eq) - 1)]))
}

cat("Computing HMM_hard...\n")
wp_hard  <- compute_wealth_path(w_hard)
RR_hard  <- wp_hard[, N_YEARS + 1] / af / salary_at_retirement
mdd_hard <- compute_mdd(wp_hard)
to_hard  <- compute_turnover(w_hard)

cat("Computing HMM_soft...\n")
wp_soft  <- compute_wealth_path(w_soft)
RR_soft  <- wp_soft[, N_YEARS + 1] / af / salary_at_retirement
mdd_soft <- compute_mdd(wp_soft)
to_soft  <- compute_turnover(w_soft)

cat(sprintf("\nHMM_hard: RR_median=%.3f, shortfall_70=%.3f, CVaR_5=%.3f, MDD=%.3f\n",
            median(RR_hard), mean(RR_hard < 0.70),
            mean(RR_hard[RR_hard <= quantile(RR_hard, 0.05)]),
            mean(mdd_hard)))
cat(sprintf("HMM_soft: RR_median=%.3f, shortfall_70=%.3f, CVaR_5=%.3f, MDD=%.3f\n\n",
            median(RR_soft), mean(RR_soft < 0.70),
            mean(RR_soft[RR_soft <= quantile(RR_soft, 0.05)]),
            mean(mdd_soft)))

build_metrics_row <- function(name, RR_vec, mdd_vec, to_vec) {
  data.frame(
    strategy = name,
    RR_mean = round(mean(RR_vec), 3),
    RR_median = round(median(RR_vec), 3),
    RR_sd = round(sd(RR_vec), 3),
    RR_p10 = round(quantile(RR_vec, 0.10), 3),
    RR_p25 = round(quantile(RR_vec, 0.25), 3),
    RR_p75 = round(quantile(RR_vec, 0.75), 3),
    RR_p90 = round(quantile(RR_vec, 0.90), 3),
    RR_iqr = round(quantile(RR_vec, 0.75) - quantile(RR_vec, 0.25), 3),
    P_shortfall_50pct = round(mean(RR_vec < 0.50), 3),
    P_shortfall_60pct = round(mean(RR_vec < 0.60), 3),
    P_shortfall_70pct = round(mean(RR_vec < 0.70), 3),
    P_shortfall_80pct = round(mean(RR_vec < 0.80), 3),
    CVaR_5pct = round(mean(RR_vec[RR_vec <= quantile(RR_vec, 0.05)]), 3),
    CVaR_10pct = round(mean(RR_vec[RR_vec <= quantile(RR_vec, 0.10)]), 3),
    RR_sharpe_like = round(mean(RR_vec) / sd(RR_vec), 3),
    MDD_mean = round(mean(mdd_vec), 3),
    MDD_median = round(median(mdd_vec), 3),
    MDD_p10 = round(quantile(mdd_vec, 0.10), 3),
    turnover_mean = round(mean(to_vec), 3),
    turnover_median = round(median(to_vec), 3),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

row_hard <- build_metrics_row("HMM_hard", RR_hard, mdd_hard, to_hard)
row_soft <- build_metrics_row("HMM_soft", RR_soft, mdd_soft, to_soft)

metrics_new <- bind_rows(metrics_old, row_hard, row_soft)
wp_new <- wp_old
wp_new[["HMM_hard"]] <- wp_hard
wp_new[["HMM_soft"]] <- wp_soft
mdd_new <- cbind(mdd_old, HMM_hard = mdd_hard, HMM_soft = mdd_soft)
to_new  <- cbind(to_old, HMM_hard = to_hard, HMM_soft = to_soft)
RR_new  <- cbind(RR_old, HMM_hard = RR_hard, HMM_soft = RR_soft)

cat("--- METRICS TABLE (FINAL, 11 STRATEGIES) ---\n")
print(metrics_new[, c("strategy", "RR_median", "RR_sd", "P_shortfall_70pct",
                      "CVaR_5pct", "MDD_mean", "turnover_mean")],
      row.names = FALSE)

saveRDS(metrics_new, "Data/metrics_full.rds")
saveRDS(wp_new, "Data/wealth_paths.rds")
saveRDS(mdd_new, "Data/max_drawdowns.rds")
saveRDS(to_new, "Data/turnovers.rds")
saveRDS(RR_new, "Data/replacement_ratios_with_hybrid_naive.rds")

cat("\n--- SAVED (11 strategies: 8 base + Hybrid + Direct_PPP + HMM_hard + HMM_soft) ---\n")
