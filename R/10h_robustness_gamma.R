# 10h_robustness_gamma.R
#
# Input:  Data/features_historical.rds
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#
# Output: Data/robustness_H_gamma.rds
#
# Sensitivity of MV strategies to the risk-aversion coefficient gamma.
# Baseline gamma = 5 (Campbell-Viceira 2002 default).
# Alternative: 3 (lower aversion, more equity) and 8 (higher aversion).
#
# Rule-based and ML strategies do not depend on gamma and are not retested.

library(dplyr)
library(tidyr)

# ---- 1. Load inputs ------------------------------------------------------
hist     <- readRDS("Data/features_historical.rds")
eq_ret   <- readRDS("Data/equity_returns.rds")
bd_ret   <- readRDS("Data/bond_returns_10y.rds")
salary   <- readRDS("Data/salary_path.rds")
af       <- readRDS("Data/annuity_factor.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- ncol(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]
EQ_CAP <- 0.70

# ---- 2. Historical monthly returns ---------------------------------------
hist <- hist %>%
  mutate(bond_ret = (lag(y10) / 100) / 12 - 9 * (y10 - lag(y10)) / 100)
ret_mat <- hist %>% select(eq_ret, bond_ret) %>%
  filter(complete.cases(.)) %>% as.matrix()

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

# ---- 3. MV solver --------------------------------------------------------
mv_optimal_weight <- function(mu, Sigma, gamma, eq_cap = 0.70) {
  w_eq_grid <- seq(0, eq_cap, by = 0.001)
  utility <- sapply(w_eq_grid, function(we) {
    w <- c(we, 1 - we)
    sum(w * mu) - (gamma / 2) * t(w) %*% Sigma %*% w
  })
  w_eq_grid[which.max(utility)]
}

# ---- 4. Wealth function --------------------------------------------------
compute_wealth_final <- function(w_eq) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 5. Loop over gamma values -------------------------------------------
gammas <- c(3, 5, 8)
all_results <- list()

cat("Computing sensitivity over risk aversion gamma...\n")

for (g in gammas) {
  cat(sprintf("  gamma = %d\n", g))
  w_p_val  <- mv_optimal_weight(mu_a,    Sig_a,    g, EQ_CAP)
  w_lw_val <- mv_optimal_weight(mu_a, Sig_lw_a,    g, EQ_CAP)
  
  w_mv_p  <- matrix(w_p_val,  N_SCEN, N_YEARS)
  w_mv_lw <- matrix(w_lw_val, N_SCEN, N_YEARS)
  
  wf_p  <- compute_wealth_final(w_mv_p)
  wf_lw <- compute_wealth_final(w_mv_lw)
  
  RR_p  <- wf_p  / af / salary_at_retirement
  RR_lw <- wf_lw / af / salary_at_retirement
  
  all_results[[paste(g, "MV_plain", sep = "_")]] <- data.frame(
    strategy = "MV_plain", gamma = g, w_equity = round(w_p_val, 3),
    RR_mean = round(mean(RR_p), 3), RR_median = round(median(RR_p), 3),
    P_shortfall_70 = round(mean(RR_p < 0.70), 3),
    stringsAsFactors = FALSE
  )
  all_results[[paste(g, "MV_LW", sep = "_")]] <- data.frame(
    strategy = "MV_LW", gamma = g, w_equity = round(w_lw_val, 3),
    RR_mean = round(mean(RR_lw), 3), RR_median = round(median(RR_lw), 3),
    P_shortfall_70 = round(mean(RR_lw < 0.70), 3),
    stringsAsFactors = FALSE
  )
}

master <- bind_rows(all_results)

cat("\n--- MV RR SENSITIVITY TO GAMMA ---\n")
print(master)

# ---- 6. Pivots -----------------------------------------------------------
cat("\n--- w_equity BY STRATEGY x GAMMA ---\n")
pivot_w <- master %>%
  select(strategy, gamma, w_equity) %>%
  pivot_wider(names_from = gamma, values_from = w_equity, names_prefix = "g=")
print(pivot_w)

cat("\n--- MEDIAN RR BY STRATEGY x GAMMA ---\n")
pivot_rr <- master %>%
  select(strategy, gamma, RR_median) %>%
  pivot_wider(names_from = gamma, values_from = RR_median, names_prefix = "g=")
print(pivot_rr)

cat("\n--- SHORTFALL_70 BY STRATEGY x GAMMA ---\n")
pivot_sf <- master %>%
  select(strategy, gamma, P_shortfall_70) %>%
  pivot_wider(names_from = gamma, values_from = P_shortfall_70, names_prefix = "g=")
print(pivot_sf)

# ---- 7. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_H_gamma.rds")
cat("\nSaved Data/robustness_H_gamma.rds\n")
