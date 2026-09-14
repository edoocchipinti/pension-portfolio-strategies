# 08_metrics.R
#
# Input:  Data/replacement_ratios_all.rds      (20000 x 8)
#         Data/dnb_equity_allocations.rds      (20000 x 42 x 5 ML/MV strategies)
#         Data/equity_returns.rds              (20000 x 100)
#         Data/bond_returns_10y.rds            (20000 x 100)
#         Data/salary_path.rds                 (20000 x 43)
#
# Output: Data/metrics_full.rds                (data.frame with all metrics)
#         Data/wealth_paths.rds                (list of 8 matrices, 20000 x 43)
#
# Dependencies: dplyr
#
# Computes extended metrics for all 8 strategies:
#   - Replacement Ratio summary (already done in 07, but re-tabulated)
#   - Shortfall probabilities at 50%, 60%, 70%, 80% RR
#   - CVaR (expected RR conditional on bottom 5% / 10%)
#   - Wealth-trajectory based: max drawdown, Sharpe-like ratio
#   - Allocation turnover (only meaningful for ML)
#
# Notes:
#   - "Sharpe-like ratio" is the cross-scenario ratio of mean RR over SD of RR,
#     since we have only ONE outcome per scenario (no time series of returns).
#     This is a measure of dispersion-adjusted performance.
#   - Max drawdown is computed on the WEALTH TRAJECTORY of each scenario
#     (peak-to-trough during accumulation), then aggregated across scenarios.

library(dplyr)

# ---- 1. Load inputs ------------------------------------------------------
RR          <- readRDS("Data/replacement_ratios_all.rds")           # 20000 x 8
allocations <- readRDS("Data/dnb_equity_allocations.rds")           # 20000 x 42 x 5
eq_ret      <- readRDS("Data/equity_returns.rds")
bd_ret      <- readRDS("Data/bond_returns_10y.rds")
salary      <- readRDS("Data/salary_path.rds")

N_SCEN  <- nrow(RR)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14

strategies <- colnames(RR)
cat("Strategies:", strategies, "\n")
cat("N_SCEN:", N_SCEN, " N_YEARS:", N_YEARS, "\n\n")

# ---- 2. Reconstruct equity-weight paths for ALL 8 strategies -------------
# Rule-based + statistical: build the equity weight matrices.
# We need these for (a) wealth trajectory, (b) turnover.

w_6040  <- matrix(0.60, nrow = N_SCEN, ncol = N_YEARS)
w_1N    <- matrix(0.50, nrow = N_SCEN, ncol = N_YEARS)
glide_path_vec <- seq(0.80, 0.30, length.out = N_YEARS)
w_glide <- matrix(rep(glide_path_vec, each = N_SCEN), nrow = N_SCEN, ncol = N_YEARS)

w_mv_plain <- allocations[, , "MV_plain"]
w_mv_lw    <- allocations[, , "MV_LW"]
w_lasso    <- allocations[, , "Lasso"]
w_ridge    <- allocations[, , "Ridge"]
w_rf       <- allocations[, , "RF"]

strategy_weights <- list(
  "60/40"    = w_6040,
  "1/N"      = w_1N,
  "Glide"    = w_glide,
  "MV_plain" = w_mv_plain,
  "MV_LW"    = w_mv_lw,
  "Lasso"    = w_lasso,
  "Ridge"    = w_ridge,
  "RF"       = w_rf
)

# ---- 3. Wealth trajectory function ---------------------------------------
# Returns the FULL wealth path (20000 x 43) so we can compute drawdown.

compute_wealth_path <- function(w_eq) {
  wealth_path <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)  # column 1 = t=0
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth_path[, t + 1] <- (wealth_path[, t] + contrib) * (1 + port_ret)
  }
  wealth_path
}

cat("Computing wealth trajectories for all 8 strategies...\n")
wealth_paths <- lapply(strategy_weights, compute_wealth_path)
names(wealth_paths) <- names(strategy_weights)

# ---- 4. Max drawdown function --------------------------------------------
# For each scenario: peak-to-trough during accumulation.
# Drawdown_t = wealth_t / max(wealth_0..t) - 1
# Max drawdown = min over t.
#
# NOTE: in a pension accumulation context with positive contributions every
# year, "drawdown" only makes sense if returns are bad enough to wipe out
# both contributions AND prior gains. We compute it as min of (wealth_t
# relative to running max) - 1.

compute_max_drawdown <- function(wealth_path) {
  # wealth_path: 20000 x 43
  # Avoid divisions by zero: at t=0, wealth=0. We start running max from t=1.
  N <- nrow(wealth_path); T <- ncol(wealth_path)
  running_max <- apply(wealth_path, 1, cummax)  # T x N (apply transposes)
  running_max <- t(running_max)                  # back to N x T
  
  # Replace 0 in running_max with NA to avoid 0/0
  running_max[running_max == 0] <- NA
  
  drawdown <- wealth_path / running_max - 1
  # Max drawdown for each scenario = min (most negative)
  apply(drawdown, 1, min, na.rm = TRUE)
}

cat("Computing max drawdown for each strategy...\n")
max_drawdowns <- sapply(wealth_paths, compute_max_drawdown)
# max_drawdowns: 20000 x 8

# ---- 5. Turnover function ------------------------------------------------
# For each scenario, sum of |Δw_t| across years.
# Higher = more rebalancing.
# Rule-based: 60/40 and 1/N have zero. Glide has small linear turnover.
# MV: zero (static). ML: variable.

compute_turnover <- function(w_eq) {
  if (ncol(w_eq) < 2) return(rep(0, nrow(w_eq)))
  delta <- abs(w_eq[, 2:ncol(w_eq)] - w_eq[, 1:(ncol(w_eq) - 1)])
  rowSums(delta)
}

cat("Computing turnover...\n")
turnovers <- sapply(strategy_weights, compute_turnover)
# turnovers: 20000 x 8

# ---- 6. Cross-scenario summary metrics -----------------------------------
# Build the master metric data.frame.

RR_TARGETS <- c(0.50, 0.60, 0.70, 0.80)
CVAR_LEVELS <- c(0.05, 0.10)

metrics <- data.frame(strategy = strategies, stringsAsFactors = FALSE)

# Replacement Ratio basics
metrics$RR_mean   <- round(colMeans(RR), 3)
metrics$RR_median <- round(apply(RR, 2, median), 3)
metrics$RR_sd     <- round(apply(RR, 2, sd), 3)
metrics$RR_p10    <- round(apply(RR, 2, quantile, 0.10), 3)
metrics$RR_p25    <- round(apply(RR, 2, quantile, 0.25), 3)
metrics$RR_p75    <- round(apply(RR, 2, quantile, 0.75), 3)
metrics$RR_p90    <- round(apply(RR, 2, quantile, 0.90), 3)
metrics$RR_iqr    <- metrics$RR_p75 - metrics$RR_p25

# Shortfall probabilities
for (tgt in RR_TARGETS) {
  col <- paste0("P_shortfall_", round(tgt * 100), "pct")
  metrics[[col]] <- round(colMeans(RR < tgt), 3)
}

# CVaR (expected RR conditional on bottom alpha %)
for (alpha in CVAR_LEVELS) {
  cvar_vals <- sapply(strategies, function(s) {
    cutoff <- quantile(RR[, s], alpha)
    mean(RR[RR[, s] <= cutoff, s])
  })
  col <- paste0("CVaR_", round(alpha * 100), "pct")
  metrics[[col]] <- round(cvar_vals, 3)
}

# Dispersion-adjusted return (Sharpe-like, cross-scenario)
# Mean RR / SD of RR. Higher = better risk-adjusted outcome dispersion.
metrics$RR_sharpe_like <- round(metrics$RR_mean / metrics$RR_sd, 3)

# Max drawdown statistics
metrics$MDD_mean    <- round(colMeans(max_drawdowns),               3)
metrics$MDD_median  <- round(apply(max_drawdowns, 2, median),       3)
metrics$MDD_p10     <- round(apply(max_drawdowns, 2, quantile, 0.10), 3)
# (p10 of drawdown = the WORST 10% drawdowns; recall MDD is negative)

# Turnover
metrics$turnover_mean   <- round(colMeans(turnovers), 3)
metrics$turnover_median <- round(apply(turnovers, 2, median), 3)

# ---- 7. Display master table ---------------------------------------------
cat("\n--- MASTER METRICS TABLE ---\n\n")
print(metrics, row.names = FALSE)

# ---- 8. Quick visual summary on key dimensions ---------------------------
cat("\n--- KEY DIMENSIONS RANKED ---\n\n")

cat("Median RR (higher = better):\n")
print(metrics %>% arrange(desc(RR_median)) %>%
        select(strategy, RR_median, RR_iqr))

cat("\nShortfall at 70% (lower = better):\n")
print(metrics %>% arrange(P_shortfall_70pct) %>%
        select(strategy, P_shortfall_70pct, P_shortfall_60pct))

cat("\nCVaR 5% (less negative = better; this is RR at the worst 5% tail):\n")
print(metrics %>% arrange(desc(CVaR_5pct)) %>%
        select(strategy, CVaR_5pct, CVaR_10pct))

cat("\nMax drawdown (less negative = better):\n")
print(metrics %>% arrange(desc(MDD_mean)) %>%
        select(strategy, MDD_mean, MDD_p10))

cat("\nTurnover (lower = less rebalancing):\n")
print(metrics %>% arrange(turnover_mean) %>%
        select(strategy, turnover_mean, turnover_median))

# ---- 9. Save -------------------------------------------------------------
saveRDS(metrics, "Data/metrics_full.rds")
saveRDS(wealth_paths, "Data/wealth_paths.rds")
saveRDS(max_drawdowns, "Data/max_drawdowns.rds")
saveRDS(turnovers, "Data/turnovers.rds")

cat("\nSaved Data/metrics_full.rds\n")
cat("Saved Data/wealth_paths.rds\n")
cat("Saved Data/max_drawdowns.rds\n")
cat("Saved Data/turnovers.rds\n")
