# 10i_robustness_contribution.R
#
# Input:  Data/dnb_equity_allocations.rds
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/mv_weights.rds
#
# Output: Data/robustness_I_contribution.rds
#
# Sensitivity of replacement ratio to contribution rate.
# Baseline contribution rate is 14% (EU Country Fiche on Pensions for NL, 2021).
# Tested values: 8%, 10%, 12%, 14%, 16%, 18%, 20%, 22%, 25%.
#
# Goal: identify the "magic contribution rate" at which the shortfall
# probability drops below 50% for each strategy.
#
# Note: this robustness test re-uses regime-based allocations from baseline
# (M0 mapping). The point is to vary the policy lever (contribution), not
# the strategy parameters.

library(dplyr)
library(tidyr)

# ---- 1. Load inputs ------------------------------------------------------
allocations  <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")
mv_weights   <- readRDS("Data/mv_weights.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
salary_at_retirement <- salary[, N_YEARS + 1]

# ---- 2. Build all 8 strategy weight matrices -----------------------------
w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
w_glide <- matrix(rep(seq(0.80, 0.30, length.out = N_YEARS),
                      each = N_SCEN), N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_weights$plain["equity"]), N_SCEN, N_YEARS)
w_mv_lw <- matrix(as.numeric(mv_weights$lw["equity"]),    N_SCEN, N_YEARS)
w_lasso <- allocations[, , "Lasso"]
w_ridge <- allocations[, , "Ridge"]
w_rf    <- allocations[, , "RF"]

strategies <- list(
  "60/40"    = w_6040,
  "1/N"      = w_1N,
  "Glide"    = w_glide,
  "MV_plain" = w_mv_p,
  "MV_LW"    = w_mv_lw,
  "Lasso"    = w_lasso,
  "Ridge"    = w_ridge,
  "RF"       = w_rf
)

# ---- 3. Wealth function with parametric contribution rate ----------------
compute_wealth_final <- function(w_eq, contrib_rate) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- contrib_rate * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 4. Sensitivity loop -------------------------------------------------
contribution_rates <- c(0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25)

cat("Computing RR for", length(contribution_rates), "contribution rates x 8 strategies...\n")
all_results <- list()

for (cr in contribution_rates) {
  cat(sprintf("  Contribution = %.0f%%\n", cr * 100))
  for (s_name in names(strategies)) {
    wf <- compute_wealth_final(strategies[[s_name]], cr)
    RR <- wf / af / salary_at_retirement
    
    all_results[[paste(cr, s_name, sep = "_")]] <- data.frame(
      strategy        = s_name,
      contrib_rate    = cr,
      RR_median       = round(median(RR), 3),
      RR_mean         = round(mean(RR), 3),
      P_shortfall_50  = round(mean(RR < 0.50), 3),
      P_shortfall_60  = round(mean(RR < 0.60), 3),
      P_shortfall_70  = round(mean(RR < 0.70), 3),
      P_shortfall_80  = round(mean(RR < 0.80), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 5. Display: Median RR pivot ----------------------------------------
cat("\n--- MEDIAN RR BY STRATEGY x CONTRIBUTION RATE ---\n")
pivot_median <- master %>%
  select(strategy, contrib_rate, RR_median) %>%
  pivot_wider(names_from = contrib_rate, values_from = RR_median,
              names_prefix = "c") %>%
  mutate(strategy = factor(strategy, levels = names(strategies))) %>%
  arrange(strategy)
print(pivot_median)

# ---- 6. Display: Shortfall at 70% (primary policy metric) ---------------
cat("\n--- SHORTFALL AT 70% BY STRATEGY x CONTRIBUTION RATE ---\n")
pivot_sf70 <- master %>%
  select(strategy, contrib_rate, P_shortfall_70) %>%
  pivot_wider(names_from = contrib_rate, values_from = P_shortfall_70,
              names_prefix = "c") %>%
  mutate(strategy = factor(strategy, levels = names(strategies))) %>%
  arrange(strategy)
print(pivot_sf70)

# ---- 7. MAGIC NUMBER: lowest contribution rate at which P_shortfall_70 ----
# ---- drops below given threshold for each strategy ---------------------
cat("\n--- MINIMUM CONTRIBUTION TO ACHIEVE THRESHOLD SHORTFALL ---\n")

magic_numbers <- function(threshold_sf70) {
  master %>%
    arrange(strategy, contrib_rate) %>%
    group_by(strategy) %>%
    summarise(
      min_contrib = {
        idx <- which(P_shortfall_70 < threshold_sf70)
        if (length(idx) > 0) contrib_rate[idx[1]] else NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      strategy = factor(strategy, levels = names(strategies)),
      threshold = paste0("sf_70_below_", threshold_sf70 * 100, "pct")
    ) %>%
    arrange(strategy)
}

cat("\nMinimum contribution rate for P(shortfall<70%) below 50%:\n")
print(magic_numbers(0.50))

cat("\nMinimum contribution rate for P(shortfall<70%) below 30%:\n")
print(magic_numbers(0.30))

cat("\nMinimum contribution rate for P(shortfall<70%) below 20%:\n")
print(magic_numbers(0.20))

# ---- 8. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_I_contribution.rds")
cat("\nSaved Data/robustness_I_contribution.rds\n")
