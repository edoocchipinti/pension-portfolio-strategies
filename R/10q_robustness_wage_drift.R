# 10q_robustness_wage_drift.R
#
# Input:  Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx
#         Data/dnb_equity_allocations.rds
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/annuity_factor.rds
#         Data/mv_weights.rds
#
# Output: Data/robustness_Q_wage_drift.rds
#
# Sensitivity of replacement ratio to the real wage drift assumption.
# Baseline: 0.4% real wage drift per year.
# Alternatives: 0.0% (no real growth) and 0.8% (more optimistic productivity).
#
# Wage drift affects:
#   - salary_path: salary[i, t] grows with scenario_inflation[i, t] + drift
#   - salary_at_retirement (final value): higher drift → larger denominator in RR
#   - contributions: higher salaries → larger contributions → larger wealth
#
# The two effects partially cancel out. The net effect on RR is non-obvious
# and depends on the relative growth of contributions vs final salary.

library(dplyr)
library(tidyr)
library(readxl)

# ---- 1. Load inputs ------------------------------------------------------
allocations  <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
af           <- readRDS("Data/annuity_factor.rds")
mv_weights   <- readRDS("Data/mv_weights.rds")

# DNB inflation (decimal, NOT scaled - we need raw decimal for wage growth)
hicp_yoy_dnb <- as.matrix(read_excel(
  "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx",
  sheet = "5_EU_Price_Inflation", col_names = FALSE))

N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14

# ---- 2. Salary path builder (parametric wage drift) ----------------------
# salary[i, 1] = 1 at age 25 (normalized)
# salary[i, t+1] = salary[i, t] * (1 + inflation[i, t] + real_drift)
# salary_path needs N_YEARS + 1 = 43 columns (year 1 to year 43 = retirement)

build_salary_path <- function(real_drift) {
  sp <- matrix(1, nrow = N_SCEN, ncol = N_YEARS + 1)
  for (t in 1:N_YEARS) {
    growth <- hicp_yoy_dnb[, t] + real_drift
    sp[, t + 1] <- sp[, t] * (1 + growth)
  }
  sp
}

# ---- 3. Strategy weight matrices (same across drift values) -------------
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
  "60/40" = w_6040, "1/N" = w_1N, "Glide" = w_glide,
  "MV_plain" = w_mv_p, "MV_LW" = w_mv_lw,
  "Lasso" = w_lasso, "Ridge" = w_ridge, "RF" = w_rf
)

# ---- 4. Wealth function (parametric salary path) ------------------------
compute_wealth_final <- function(w_eq, salary) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 5. Sensitivity loop -------------------------------------------------
drifts <- c(0.000, 0.004, 0.008)
labels <- c("drift_0pct", "drift_0.4pct_baseline", "drift_0.8pct")
all_results <- list()

cat("Computing sensitivity over real wage drift...\n")

for (i in seq_along(drifts)) {
  d   <- drifts[i]
  lbl <- labels[i]
  cat(sprintf("  Drift = %.1f%% (%s)\n", 100 * d, lbl))
  
  salary <- build_salary_path(d)
  salary_at_retirement <- salary[, N_YEARS + 1]
  
  for (s_name in names(strategies)) {
    wf <- compute_wealth_final(strategies[[s_name]], salary)
    RR <- wf / af / salary_at_retirement
    all_results[[paste(lbl, s_name, sep = "_")]] <- data.frame(
      strategy = s_name, drift = lbl, drift_val = d,
      RR_mean = round(mean(RR), 3),
      RR_median = round(median(RR), 3),
      RR_sd = round(sd(RR), 3),
      P_shortfall_70 = round(mean(RR < 0.70), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 6. Display pivots ---------------------------------------------------
strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF")

cat("\n--- MEDIAN RR BY STRATEGY x WAGE DRIFT ---\n")
pivot_median <- master %>%
  select(strategy, drift, RR_median) %>%
  pivot_wider(names_from = drift, values_from = RR_median) %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_median)

cat("\n--- SHORTFALL_70 BY STRATEGY x WAGE DRIFT ---\n")
pivot_sf <- master %>%
  select(strategy, drift, P_shortfall_70) %>%
  pivot_wider(names_from = drift, values_from = P_shortfall_70) %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_sf)

# ---- 7. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_Q_wage_drift.rds")
cat("\nSaved Data/robustness_Q_wage_drift.rds\n")
