# 10d_robustness_bond_duration.R
#
# Input:  Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx
#         Data/yield_10y.rds
#         Data/dnb_equity_allocations.rds
#         Data/equity_returns.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/mv_weights.rds
#
# Output: Data/robustness_D_bond_duration.rds
#
# Sensitivity of replacement ratio to the duration assumption used to
# construct bond returns. The default rolling 10-year bond is approximated
# with effective duration ~9 (D=9), matching the standard pension-fund
# accounting convention. Alternative durations D=8 and D=10 test the
# robustness of the rule-based, MV, and ML strategy rankings.
#
# Method:
#   bond_return_t = -D * (y10_{t+1} - y10_{t}) + y10_t
#                  (first-order duration approximation, in decimal)
#
# Reuses ML class predictions and MV weights (independent of bond duration).

library(dplyr)
library(tidyr)

# ---- 1. Load inputs ------------------------------------------------------
y10          <- as.matrix(readRDS("Data/yield_10y.rds"))     # decimal, 20000 x 101
eq_ret       <- readRDS("Data/equity_returns.rds")
allocations  <- readRDS("Data/dnb_equity_allocations.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")
mv_weights   <- readRDS("Data/mv_weights.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

# ---- 2. Bond return constructor ------------------------------------------
# Duration approximation: r_bond_t ≈ y_t - D * Δy_t
# where Δy_t = y_{t+1} - y_t (yield realized at end of year t+1 vs start).
# Uses yields at t (column t in DNB indexing) and t+1.

compute_bond_returns <- function(D) {
  # y10 has 101 columns (t=0 to t=100). Bond returns for year t use
  # y10[, t] and y10[, t+1]. We compute the first N_YEARS+1 columns to
  # match equity_returns' length (100), but only use first N_YEARS.
  br <- matrix(NA, nrow = N_SCEN, ncol = ncol(y10) - 1)
  for (t in 1:(ncol(y10) - 1)) {
    y_now  <- y10[, t]
    y_next <- y10[, t + 1]
    br[, t] <- y_now - D * (y_next - y_now)
  }
  br
}

# ---- 3. Strategy weights matrices (same for all bond durations) ----------
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

# ---- 4. Wealth function (now bond returns are parametric) ---------------
compute_wealth_final <- function(w_eq, bd_ret) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 5. Loop over durations ---------------------------------------------
durations <- c(8, 9, 10)
all_results <- list()

cat("Computing sensitivity over bond duration...\n")
for (D in durations) {
  cat(sprintf("  Duration D=%d\n", D))
  bd_ret <- compute_bond_returns(D)
  
  for (s_name in names(strategies)) {
    wf <- compute_wealth_final(strategies[[s_name]], bd_ret)
    RR <- wf / af / salary_at_retirement
    all_results[[paste(D, s_name, sep = "_")]] <- data.frame(
      strategy        = s_name,
      duration        = D,
      RR_mean         = round(mean(RR), 3),
      RR_median       = round(median(RR), 3),
      RR_sd           = round(sd(RR), 3),
      P_shortfall_70  = round(mean(RR < 0.70), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 6. Display pivots ---------------------------------------------------
strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF")

cat("\n--- MEDIAN RR BY STRATEGY x DURATION ---\n")
pivot_median <- master %>%
  select(strategy, duration, RR_median) %>%
  pivot_wider(names_from = duration, values_from = RR_median,
              names_prefix = "D=") %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_median)

cat("\n--- SHORTFALL_70 BY STRATEGY x DURATION ---\n")
pivot_sf <- master %>%
  select(strategy, duration, P_shortfall_70) %>%
  pivot_wider(names_from = duration, values_from = P_shortfall_70,
              names_prefix = "D=") %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_sf)

# ---- 7. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_D_bond_duration.rds")
cat("\nSaved Data/robustness_D_bond_duration.rds\n")
