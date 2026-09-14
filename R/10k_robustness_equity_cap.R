# 10k_robustness_equity_cap.R
#
# Input:  Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/features_historical.rds       (for MV recalibration)
#         Data/ml_models.rds                  (re-used for ML predictions)
#         Data/feature_cols.rds
#         Data/yield_*.rds
#         Data/cp2022-...xlsx
#
# Output: Data/robustness_K_equity_cap.rds
#
# Sensitivity of replacement ratio to the IORP II equity cap.
# Baseline cap is 70% (IORP II Article 19).
# Tested caps: 60%, 70%, 80%, 90%, 100%.
#
# Effects:
#   - 60/40, 1/N: unaffected (their equity weights are 60%, 50%)
#   - Glide-path: starting point capped at min(0.80, cap); end point 0.30 unchanged
#   - MV: re-optimised with the new cap (binding constraint)
#   - ML: regime->equity mapping CAPPED at new cap. So baseline 70/55/45/30
#         becomes min(cap, 70)/min(cap,55)/.../min(cap,30).
#
# Note: this is a regulatory "what-if". IORP II is a hard constraint in
# practice. The robustness test shows how results would change if the cap
# were relaxed, NOT a recommendation to relax it.

library(dplyr)
library(tidyr)
library(glmnet)
library(ranger)
library(readxl)

# ---- 1. Load inputs ------------------------------------------------------
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")
hist         <- readRDS("Data/features_historical.rds")
allocations  <- readRDS("Data/dnb_equity_allocations.rds")  # for ML class labels

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]
GAMMA <- 5

# ---- 2. Reload ML class predictions from baseline allocations ----------
# The baseline allocations file already contains M0 mapping (70/55/45/30).
# We need to recover the underlying CLASS predictions to re-map them with
# different caps. Since we don't store classes directly, we use a reverse
# lookup based on the unique equity weights.
#
# Baseline map M0: 70=Normal, 55=Inflationary, 45=Late, 30=Stagflation
baseline_map <- c(Normal_expansion = 0.70, Inflationary_exp = 0.55,
                  Late_cycle = 0.45, Stagflation_risk = 0.30)
reverse_lookup <- function(equity_weight) {
  names(baseline_map)[match(round(equity_weight, 2), round(baseline_map, 2))]
}

ml_classes <- list(
  Lasso = matrix(reverse_lookup(allocations[, , "Lasso"]),
                 N_SCEN, N_YEARS),
  Ridge = matrix(reverse_lookup(allocations[, , "Ridge"]),
                 N_SCEN, N_YEARS),
  RF    = matrix(reverse_lookup(allocations[, , "RF"]),
                 N_SCEN, N_YEARS)
)

# ---- 3. MV recalibration with arbitrary cap -----------------------------
ret_mat <- as.matrix(hist[, c("eq_ret")])
hist_bond_proxy <- (lag(hist$y10) / 100) / 12 - 9 * (hist$y10 - lag(hist$y10)) / 100
ret_mat <- cbind(eq_ret = hist$eq_ret, bond_ret = hist_bond_proxy)
ret_mat <- ret_mat[complete.cases(ret_mat), ]

mu_a  <- colMeans(ret_mat) * 12
Sig_a <- cov(ret_mat) * 12

ledoit_wolf_shrinkage <- function(R) {
  T <- nrow(R); N <- ncol(R)
  Rc <- scale(R, center = TRUE, scale = FALSE)
  S <- cov(R)
  F <- diag(mean(diag(S)), N)
  pi_mat <- (t(Rc^2) %*% (Rc^2)) / T - S^2
  pi_hat <- sum(pi_mat); rho_hat <- sum(diag(pi_mat))
  gamma_hat <- sum((S - F)^2)
  delta <- max(0, min(1, (pi_hat - rho_hat) / gamma_hat / T))
  delta * F + (1 - delta) * S
}
Sig_lw_a <- ledoit_wolf_shrinkage(ret_mat) * 12

mv_optimal_weights <- function(mu, Sigma, gamma = 5, eq_cap = 0.70) {
  w_eq_grid <- seq(0, eq_cap, by = 0.001)
  utility <- sapply(w_eq_grid, function(we) {
    w <- c(we, 1 - we)
    sum(w * mu) - (gamma / 2) * t(w) %*% Sigma %*% w
  })
  w_eq_grid[which.max(utility)]
}

# ---- 4. Wealth function -------------------------------------------------
compute_wealth_final <- function(w_eq) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 5. Sensitivity loop ------------------------------------------------
caps <- c(0.60, 0.70, 0.80, 0.90, 1.00)
all_results <- list()

cat("Computing sensitivity over equity caps...\n")

for (cap in caps) {
  cat(sprintf("  Cap = %.0f%%\n", cap * 100))
  
  # 60/40, 1/N unaffected
  w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
  w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
  
  # Glide: start = min(0.80, cap), end = 0.30
  glide_start <- min(0.80, cap)
  w_glide <- matrix(rep(seq(glide_start, 0.30, length.out = N_YEARS),
                        each = N_SCEN), N_SCEN, N_YEARS)
  
  # MV: re-optimise with new cap
  w_mv_plain_val <- mv_optimal_weights(mu_a, Sig_a, GAMMA, cap)
  w_mv_lw_val    <- mv_optimal_weights(mu_a, Sig_lw_a, GAMMA, cap)
  w_mv_p  <- matrix(w_mv_plain_val, N_SCEN, N_YEARS)
  w_mv_lw <- matrix(w_mv_lw_val,    N_SCEN, N_YEARS)
  
  # ML: capped allocation map (cap baseline at new cap)
  capped_map <- pmin(baseline_map, cap)
  w_lasso <- matrix(capped_map[ml_classes$Lasso], N_SCEN, N_YEARS)
  w_ridge <- matrix(capped_map[ml_classes$Ridge], N_SCEN, N_YEARS)
  w_rf    <- matrix(capped_map[ml_classes$RF],    N_SCEN, N_YEARS)
  
  strategies <- list(
    "60/40" = w_6040, "1/N" = w_1N, "Glide" = w_glide,
    "MV_plain" = w_mv_p, "MV_LW" = w_mv_lw,
    "Lasso" = w_lasso, "Ridge" = w_ridge, "RF" = w_rf
  )
  
  for (s_name in names(strategies)) {
    wf <- compute_wealth_final(strategies[[s_name]])
    RR <- wf / af / salary_at_retirement
    all_results[[paste(cap, s_name, sep = "_")]] <- data.frame(
      strategy        = s_name,
      cap             = cap,
      mean_eq_weight  = round(mean(strategies[[s_name]]), 3),
      RR_median       = round(median(RR), 3),
      RR_mean         = round(mean(RR), 3),
      P_shortfall_70  = round(mean(RR < 0.70), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 6. Display pivots ---------------------------------------------------
strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF")

cat("\n--- MEDIAN RR BY STRATEGY x EQUITY CAP ---\n")
pivot_median <- master %>%
  select(strategy, cap, RR_median) %>%
  pivot_wider(names_from = cap, values_from = RR_median, names_prefix = "cap_") %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_median)

cat("\n--- SHORTFALL_70 BY STRATEGY x EQUITY CAP ---\n")
pivot_sf <- master %>%
  select(strategy, cap, P_shortfall_70) %>%
  pivot_wider(names_from = cap, values_from = P_shortfall_70, names_prefix = "cap_") %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_sf)

cat("\n--- MEAN EQUITY WEIGHT BY STRATEGY x CAP ---\n")
pivot_w <- master %>%
  select(strategy, cap, mean_eq_weight) %>%
  pivot_wider(names_from = cap, values_from = mean_eq_weight, names_prefix = "cap_") %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_w)

# ---- 7. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_K_equity_cap.rds")
cat("\nSaved Data/robustness_K_equity_cap.rds\n")

# 10L_robustness_target.R
#
# Input:  Data/replacement_ratios_all.rds   (20000 x 8)
#
# Output: Data/robustness_L_target.rds
#
# Sensitivity of the shortfall metric to alternative replacement ratio
# targets. Baseline WTP target is 70%. Tested values: 50%, 55%, 60%, 65%,
# 70%, 75%, 80%.
#
# This robustness has zero re-computation cost: we just re-tabulate
# shortfall against different thresholds on the existing RR matrix.
#
# Goal: show how the ranking and absolute shortfall changes if the
# regulatory or fund-level target were less ambitious (60%, similar to
# Germany or UK) or more ambitious (80%).

library(dplyr)
library(tidyr)

# ---- 1. Load --------------------------------------------------------------
RR <- readRDS("Data/replacement_ratios_all.rds")

strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF")
targets <- c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80)

# ---- 2. Compute shortfall against each target ----------------------------
all_results <- list()

for (tgt in targets) {
  for (s in strategy_order) {
    all_results[[paste(tgt, s, sep = "_")]] <- data.frame(
      strategy        = s,
      target          = tgt,
      shortfall       = round(mean(RR[, s] < tgt), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 3. Display: shortfall pivot -----------------------------------------
cat("--- SHORTFALL PROBABILITY BY STRATEGY x TARGET ---\n\n")
pivot <- master %>%
  pivot_wider(names_from = target, values_from = shortfall,
              names_prefix = "tgt_") %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot)

# ---- 4. Identify ranking changes ----------------------------------------
cat("\n--- RANKING (best = lowest shortfall) BY TARGET ---\n\n")
ranking <- master %>%
  group_by(target) %>%
  arrange(target, shortfall) %>%
  mutate(rank = row_number()) %>%
  select(target, rank, strategy, shortfall) %>%
  pivot_wider(names_from = rank, values_from = strategy,
              names_prefix = "rank") %>%
  ungroup()
print(ranking)

# ---- 5. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_L_target.rds")
cat("\nSaved Data/robustness_L_target.rds\n")

# 10L_robustness_target.R
#
# Input:  Data/replacement_ratios_all.rds  (RR matrix 20000 x 8)
#
# Output: Data/robustness_L_target.rds
#
# Sensitivity of shortfall probability to the pension target choice.
# Baseline target: 70% (WTP convention in the Netherlands).
# Alternatives: 50% (UK style), 60% (Germany style), 80% (more demanding).

library(dplyr)
library(tidyr)

# ---- 1. Load inputs ------------------------------------------------------
RR <- readRDS("Data/replacement_ratios_all.rds")

cat("Strategies in RR matrix:\n")
print(colnames(RR))
cat(sprintf("\nMatrix dim: %d scenarios x %d strategies\n\n",
            nrow(RR), ncol(RR)))

# ---- 2. Compute shortfall for each target -------------------------------
targets <- c(0.50, 0.60, 0.70, 0.80)
target_labels <- c("T_50pct", "T_60pct_DE", "T_70pct_baseline_WTP",
                   "T_80pct_demanding")
all_results <- list()

for (i in seq_along(targets)) {
  tg <- targets[i]
  lbl <- target_labels[i]
  cat(sprintf("Target = %.0f%% (%s)\n", 100 * tg, lbl))
  
  for (s in colnames(RR)) {
    RR_vec <- RR[, s]
    all_results[[paste(lbl, s, sep = "_")]] <- data.frame(
      strategy        = s,
      target          = lbl,
      target_val      = tg,
      RR_median       = round(median(RR_vec), 3),
      RR_mean         = round(mean(RR_vec), 3),
      P_shortfall     = round(mean(RR_vec < tg), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 3. Display pivots ---------------------------------------------------
strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF")

cat("\n--- SHORTFALL BY STRATEGY x TARGET ---\n")
pivot_sf <- master %>%
  select(strategy, target, P_shortfall) %>%
  pivot_wider(names_from = target, values_from = P_shortfall) %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_sf)

cat("\n--- BEST STRATEGY (LOWEST SHORTFALL) BY TARGET ---\n")
best_by_target <- master %>%
  group_by(target) %>%
  filter(P_shortfall == min(P_shortfall)) %>%
  select(target, target_val, strategy, P_shortfall)
print(best_by_target)

# ---- 4. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_L_target.rds")
cat("\nSaved Data/robustness_L_target.rds\n")
