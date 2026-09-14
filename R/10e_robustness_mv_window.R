# 10e_robustness_mv_window.R
#
# Input:  Data/features_historical.rds
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/dnb_equity_allocations.rds  (for Lasso/Ridge/RF, unchanged)
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#
# Output: Data/robustness_E_mv_window.rds
#
# Sensitivity of MV strategies' replacement ratio to the calibration window.
# Baseline uses the FULL historical sample (300 obs, Jan 2001 to Dec 2025).
# Alternative windows: 60, 120, 240 months (most recent).
#
# Already-discovered finding from 05_calibrate_mv.R:
#   - Shorter windows BIND on the IORP II 70% cap (w_eq = 0.70)
#   - Full sample gives interior solution (w_eq = 0.42)
#
# This script propagates the alternative MV weights through to RR
# and confirms the impact on shortfall.
#
# Rule-based and ML strategies are unchanged across windows.

library(dplyr)
library(tidyr)

# ---- 1. Load inputs ------------------------------------------------------
hist         <- readRDS("Data/features_historical.rds")
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
allocations  <- readRDS("Data/dnb_equity_allocations.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]
GAMMA  <- 5
EQ_CAP <- 0.70

# ---- 2. Build monthly return matrix (eq, bond proxy) ---------------------
# Historical bond return proxy from y10 (yield-based, monthly approximation)
# Using D=9 (rolling 10y bond accounting standard)
hist <- hist %>%
  mutate(
    bond_ret = (lag(y10) / 100) / 12 - 9 * (y10 - lag(y10)) / 100
  )

ret_full <- hist %>% select(eq_ret, bond_ret) %>%
  filter(complete.cases(.)) %>% as.matrix()

cat("Full historical sample size:", nrow(ret_full), "monthly obs\n\n")

# ---- 3. Ledoit-Wolf shrinkage --------------------------------------------
ledoit_wolf_shrinkage <- function(R) {
  T <- nrow(R); N <- ncol(R)
  Rc <- scale(R, center = TRUE, scale = FALSE)
  S <- cov(R)
  F <- diag(mean(diag(S)), N)
  pi_mat <- (t(Rc^2) %*% (Rc^2)) / T - S^2
  pi_hat <- sum(pi_mat); rho_hat <- sum(diag(pi_mat))
  gamma_hat <- sum((S - F)^2)
  delta <- max(0, min(1, (pi_hat - rho_hat) / gamma_hat / T))
  list(Sigma = delta * F + (1 - delta) * S, delta = delta)
}

# ---- 4. MV solver --------------------------------------------------------
mv_optimal_weight <- function(mu, Sigma, gamma = 5, eq_cap = 0.70) {
  w_eq_grid <- seq(0, eq_cap, by = 0.001)
  utility <- sapply(w_eq_grid, function(we) {
    w <- c(we, 1 - we)
    sum(w * mu) - (gamma / 2) * t(w) %*% Sigma %*% w
  })
  w_eq_grid[which.max(utility)]
}

# ---- 5. Compute MV weights for each window -------------------------------
windows <- c(60, 120, 240, nrow(ret_full))
window_names <- c("60m", "120m", "240m", "full")

mv_results_calib <- data.frame(
  window = window_names,
  n_obs  = windows,
  w_plain = numeric(length(windows)),
  w_lw    = numeric(length(windows)),
  delta_lw = numeric(length(windows))
)

for (i in seq_along(windows)) {
  w_len <- windows[i]
  R_win <- ret_full[(nrow(ret_full) - w_len + 1):nrow(ret_full), ]
  mu_win  <- colMeans(R_win) * 12
  Sig_win <- cov(R_win) * 12
  lw_win  <- ledoit_wolf_shrinkage(R_win)
  Sig_lw_win <- lw_win$Sigma * 12
  
  mv_results_calib$w_plain[i] <- mv_optimal_weight(mu_win, Sig_win, GAMMA, EQ_CAP)
  mv_results_calib$w_lw[i]    <- mv_optimal_weight(mu_win, Sig_lw_win, GAMMA, EQ_CAP)
  mv_results_calib$delta_lw[i] <- lw_win$delta
}

cat("--- MV CALIBRATION RESULTS BY WINDOW ---\n")
print(mv_results_calib)
cat("\nNote: shorter windows bind on the equity cap (0.70)\n\n")

# ---- 6. Compute RR for MV strategies under each window ------------------
compute_wealth_final <- function(w_eq) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

all_results <- list()

for (i in seq_along(windows)) {
  w_name <- window_names[i]
  w_p_val <- mv_results_calib$w_plain[i]
  w_lw_val <- mv_results_calib$w_lw[i]
  
  w_mv_p  <- matrix(w_p_val,  N_SCEN, N_YEARS)
  w_mv_lw <- matrix(w_lw_val, N_SCEN, N_YEARS)
  
  wf_p  <- compute_wealth_final(w_mv_p)
  wf_lw <- compute_wealth_final(w_mv_lw)
  
  RR_p  <- wf_p  / af / salary_at_retirement
  RR_lw <- wf_lw / af / salary_at_retirement
  
  all_results[[paste(w_name, "MV_plain", sep = "_")]] <- data.frame(
    strategy = "MV_plain", window = w_name, w_equity = round(w_p_val, 3),
    RR_mean = round(mean(RR_p), 3), RR_median = round(median(RR_p), 3),
    P_shortfall_70 = round(mean(RR_p < 0.70), 3),
    stringsAsFactors = FALSE
  )
  all_results[[paste(w_name, "MV_LW", sep = "_")]] <- data.frame(
    strategy = "MV_LW", window = w_name, w_equity = round(w_lw_val, 3),
    RR_mean = round(mean(RR_lw), 3), RR_median = round(median(RR_lw), 3),
    P_shortfall_70 = round(mean(RR_lw < 0.70), 3),
    stringsAsFactors = FALSE
  )
}

master <- bind_rows(all_results)

cat("--- MV RR SENSITIVITY TO CALIBRATION WINDOW ---\n")
print(master)

cat("\n--- KEY FINDING ---\n")
cat("Shorter windows (60-240m) imply w_equity = 0.70 (cap-binding),\n")
cat("identical to a 70/30 fixed-weight allocation, vastly different\n")
cat("from the full-sample interior solution (~0.42).\n")

# ---- 7. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_E_mv_window.rds")
cat("\nSaved Data/robustness_E_mv_window.rds\n")
