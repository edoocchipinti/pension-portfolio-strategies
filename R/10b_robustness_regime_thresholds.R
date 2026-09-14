# 10b_robustness_regime_thresholds.R
#
# Input:  Data/features_historical.rds
#         Data/feature_cols.rds
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/yield_*.rds
#         Data/cp2022-...xlsx
#         Data/mv_weights.rds
#
# Output: Data/robustness_B_regime_thresholds.rds
#
# Sensitivity of ML strategies to the regime taxonomy thresholds.
# Baseline (M0):
#   slope_thr = 1.0% (sample median)
#   infl_thr  = 2.0% (ECB target)
#   dd_thr    = 15% (Faber 2007)
#
# Alternatives:
#   - slope Q1 (0.4%) or Q3 (1.7%)  -- baseline uses median (1.0%)
#   - inflation 1.5% or 2.5%        -- baseline 2.0%
#
# For each alternative combination, re-label regimes, retrain Lasso/Ridge/RF,
# apply to DNB, compute RR.
#
# Combinations tested (4 alternatives + baseline = 5 cells):
#   slope_thr × infl_thr in {(0.4, 2.0), (1.7, 2.0), (1.0, 1.5), (1.0, 2.5), baseline}
#
# Note: drawdown threshold is the SAME (15%) in all cells. Tested in 10c.
#
# Total runtime: ~5 minutes (4 retrainings + 4 DNB inferences).

library(dplyr)
library(tidyr)
library(glmnet)
library(ranger)
library(readxl)

set.seed(42)

# ---- 1. Load inputs ------------------------------------------------------
hist         <- readRDS("Data/features_historical.rds")
feature_cols <- readRDS("Data/feature_cols.rds")
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")
mv_weights   <- readRDS("Data/mv_weights.rds")

y1  <- as.matrix(readRDS("Data/yield_1y.rds"))  * 100
y5  <- as.matrix(readRDS("Data/yield_5y.rds"))  * 100
y10 <- as.matrix(readRDS("Data/yield_10y.rds")) * 100
y20 <- as.matrix(readRDS("Data/yield_20y.rds")) * 100
hicp_yoy_dnb <- as.matrix(read_excel(
  "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx",
  sheet = "5_EU_Price_Inflation", col_names = FALSE)) * 100

hist_last <- tail(hist, 1)
warm <- list(
  eq_ret_l1 = hist_last$eq_ret, eq_ret_l3 = hist_last$eq_ret_l3,
  eq_ret_l12 = hist_last$eq_ret_l12,
  hicp_yoy_l3 = hist_last$hicp_yoy_l3, hicp_yoy_l12 = hist_last$hicp_yoy_l12,
  slope_l3 = hist_last$slope_l3, slope_l12 = hist_last$slope_l12,
  y10_l12 = hist_last$y10_l12
)

N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]
DD_THR <- -0.15

alloc_map <- c(Normal_expansion = 0.70, Inflationary_exp = 0.55,
               Late_cycle = 0.45, Stagflation_risk = 0.30)

# ---- 2. Regime relabel function ------------------------------------------
relabel_regimes <- function(df, slope_thr, infl_thr, dd_thr = -0.15) {
  df %>%
    mutate(
      regime = case_when(
        eq_drawdown_12m <= dd_thr & hicp_yoy >= infl_thr ~ "Stagflation_risk",
        eq_drawdown_12m <= dd_thr & hicp_yoy <  infl_thr ~ "Late_cycle",
        slope_10_1 >= slope_thr & hicp_yoy <  infl_thr   ~ "Normal_expansion",
        slope_10_1 >= slope_thr & hicp_yoy >= infl_thr   ~ "Inflationary_exp",
        slope_10_1 <  slope_thr & hicp_yoy <  infl_thr   ~ "Late_cycle",
        slope_10_1 <  slope_thr & hicp_yoy >= infl_thr   ~ "Stagflation_risk",
        TRUE ~ NA_character_
      )
    )
}

# ---- 3. ML train function ------------------------------------------------
train_ml_models <- function(df_labeled) {
  X <- as.matrix(df_labeled[, feature_cols])
  y <- factor(df_labeled$regime,
              levels = c("Normal_expansion", "Inflationary_exp",
                         "Late_cycle", "Stagflation_risk"))
  
  cat("    Training Lasso...\n")
  fit_lasso <- cv.glmnet(X, y, family = "multinomial",
                         alpha = 1, type.measure = "class")
  cat("    Training Ridge...\n")
  fit_ridge <- cv.glmnet(X, y, family = "multinomial",
                         alpha = 0, type.measure = "class")
  cat("    Training RF...\n")
  fit_rf <- ranger(x = X, y = y, num.trees = 500, probability = TRUE,
                   seed = 42)
  
  list(lasso = fit_lasso, ridge = fit_ridge, rf = fit_rf)
}

# ---- 4. DNB feature builder ----------------------------------------------
build_features_year <- function(t_career) {
  t_yc <- t_career + 1; t_yp <- t_career
  t_rc <- t_career; t_rp <- max(t_career - 1, 1)
  t_ic <- t_career; t_ip <- max(t_career - 1, 1)
  
  eq_ret_t   <- eq_ret[, t_rc]
  hicp_yoy_t <- hicp_yoy_dnb[, t_ic]
  y1_t <- y1[, t_yc]; y5_t <- y5[, t_yc]
  y10_t <- y10[, t_yc]; y20_t <- y20[, t_yc]
  
  if (t_career == 1) { y1_prev <- y1[, 1]; y10_prev <- y10[, 1] }
  else { y1_prev <- y1[, t_yp]; y10_prev <- y10[, t_yp] }
  dy1_t  <- 100 * (y1_t  - y1_prev)
  dy10_t <- 100 * (y10_t - y10_prev)
  
  slope_10_1_t <- y10_t - y1_t
  slope_20_5_t <- y20_t - y5_t
  curvature_t  <- 2 * y5_t - y1_t - y10_t
  real_yield10_t <- y10_t - hicp_yoy_t
  eq_drawdown_12m_t <- pmin(0, eq_ret_t)
  
  if (t_career == 1) {
    eq_ret_l1_t    <- rep(warm$eq_ret_l1,    N_SCEN)
    eq_ret_l3_t    <- rep(warm$eq_ret_l3,    N_SCEN)
    eq_ret_l12_t   <- rep(warm$eq_ret_l12,   N_SCEN)
    hicp_yoy_l3_t  <- rep(warm$hicp_yoy_l3,  N_SCEN)
    hicp_yoy_l12_t <- rep(warm$hicp_yoy_l12, N_SCEN)
    slope_l3_t     <- rep(warm$slope_l3,     N_SCEN)
    slope_l12_t    <- rep(warm$slope_l12,    N_SCEN)
    y10_l12_t      <- rep(warm$y10_l12,      N_SCEN)
  } else {
    eq_ret_l1_t   <- eq_ret[, t_rp]; eq_ret_l3_t <- eq_ret_l1_t
    eq_ret_l12_t  <- eq_ret[, t_rp]
    hicp_yoy_l3_t <- hicp_yoy_t
    hicp_yoy_l12_t <- hicp_yoy_dnb[, t_ip]
    y1p <- y1[, t_yp]; y10p <- y10[, t_yp]
    sp  <- y10p - y1p
    slope_l3_t <- sp; slope_l12_t <- sp; y10_l12_t <- y10p
  }
  slope_x_infl_t <- slope_10_1_t * hicp_yoy_t
  realy_x_eq_t   <- real_yield10_t * eq_ret_l1_t
  
  X <- cbind(eq_ret = eq_ret_t, hicp_yoy = hicp_yoy_t,
             y1 = y1_t, y5 = y5_t, y10 = y10_t, y20 = y20_t,
             dy1 = dy1_t, dy10 = dy10_t,
             slope_10_1 = slope_10_1_t, slope_20_5 = slope_20_5_t,
             curvature = curvature_t, real_yield10 = real_yield10_t,
             eq_drawdown_12m = eq_drawdown_12m_t,
             eq_ret_l1 = eq_ret_l1_t, eq_ret_l3 = eq_ret_l3_t,
             eq_ret_l12 = eq_ret_l12_t,
             hicp_yoy_l3 = hicp_yoy_l3_t, hicp_yoy_l12 = hicp_yoy_l12_t,
             slope_l3 = slope_l3_t, slope_l12 = slope_l12_t,
             y10_l12 = y10_l12_t,
             slope_x_infl = slope_x_infl_t, realy_x_eq = realy_x_eq_t)
  X
}

# ---- 5. DNB inference function -------------------------------------------
apply_to_dnb <- function(models) {
  allocs <- array(NA_real_, dim = c(N_SCEN, N_YEARS, 3),
                  dimnames = list(NULL, NULL, c("Lasso", "Ridge", "RF")))
  for (t in 1:N_YEARS) {
    X_t <- build_features_year(t)
    pl <- predict(models$lasso, newx = X_t, s = "lambda.min", type = "response")[, , 1]
    allocs[, t, "Lasso"] <- alloc_map[colnames(pl)[apply(pl, 1, which.max)]]
    pr <- predict(models$ridge, newx = X_t, s = "lambda.min", type = "response")[, , 1]
    allocs[, t, "Ridge"] <- alloc_map[colnames(pr)[apply(pr, 1, which.max)]]
    pf <- predict(models$rf, data = as.data.frame(X_t))$predictions
    allocs[, t, "RF"] <- alloc_map[colnames(pf)[apply(pf, 1, which.max)]]
  }
  allocs
}

# ---- 6. Wealth and RR ----------------------------------------------------
compute_wealth_final <- function(w_eq) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# ---- 7. Compute RR for the 5 invariant strategies (once) ----------------
w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
w_glide <- matrix(rep(seq(0.80, 0.30, length.out = N_YEARS),
                      each = N_SCEN), N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_weights$plain["equity"]), N_SCEN, N_YEARS)
w_mv_lw <- matrix(as.numeric(mv_weights$lw["equity"]),    N_SCEN, N_YEARS)

invariant_RR <- function() {
  list(
    `60/40`    = compute_wealth_final(w_6040)  / af / salary_at_retirement,
    `1/N`      = compute_wealth_final(w_1N)    / af / salary_at_retirement,
    Glide      = compute_wealth_final(w_glide) / af / salary_at_retirement,
    MV_plain   = compute_wealth_final(w_mv_p)  / af / salary_at_retirement,
    MV_LW      = compute_wealth_final(w_mv_lw) / af / salary_at_retirement
  )
}

cat("Computing invariant strategy RR (rule-based + MV)...\n")
inv_RR <- invariant_RR()

# ---- 8. Main loop over threshold combinations ----------------------------
combinations <- list(
  baseline    = list(slope_thr = 1.0, infl_thr = 2.0, label = "M0_baseline"),
  slope_q1    = list(slope_thr = 0.4, infl_thr = 2.0, label = "Slope_Q1"),
  slope_q3    = list(slope_thr = 1.7, infl_thr = 2.0, label = "Slope_Q3"),
  infl_low    = list(slope_thr = 1.0, infl_thr = 1.5, label = "Infl_1.5pct"),
  infl_high   = list(slope_thr = 1.0, infl_thr = 2.5, label = "Infl_2.5pct")
)

all_results <- list()

for (cmb_name in names(combinations)) {
  cmb <- combinations[[cmb_name]]
  cat(sprintf("\n=== Combination: %s (slope=%.1f, infl=%.1f) ===\n",
              cmb$label, cmb$slope_thr, cmb$infl_thr))
  
  # Relabel
  df_lab <- relabel_regimes(hist, cmb$slope_thr, cmb$infl_thr, DD_THR)
  cat("  Regime distribution:\n")
  print(table(df_lab$regime))
  
  # Train ML
  models_cmb <- train_ml_models(df_lab)
  
  # Apply to DNB
  cat("  Applying to DNB...\n")
  allocs <- apply_to_dnb(models_cmb)
  
  # Compute RR for ML and combine with invariant
  RR_lasso <- compute_wealth_final(allocs[, , "Lasso"]) / af / salary_at_retirement
  RR_ridge <- compute_wealth_final(allocs[, , "Ridge"]) / af / salary_at_retirement
  RR_rf    <- compute_wealth_final(allocs[, , "RF"])    / af / salary_at_retirement
  
  RR_all <- list(
    `60/40` = inv_RR$`60/40`, `1/N` = inv_RR$`1/N`, Glide = inv_RR$Glide,
    MV_plain = inv_RR$MV_plain, MV_LW = inv_RR$MV_LW,
    Lasso = RR_lasso, Ridge = RR_ridge, RF = RR_rf
  )
  
  for (s_name in names(RR_all)) {
    rr_vec <- RR_all[[s_name]]
    all_results[[paste(cmb_name, s_name, sep = "_")]] <- data.frame(
      strategy = s_name, combination = cmb$label,
      slope_thr = cmb$slope_thr, infl_thr = cmb$infl_thr,
      RR_median = round(median(rr_vec), 3),
      RR_mean = round(mean(rr_vec), 3),
      P_shortfall_70 = round(mean(rr_vec < 0.70), 3),
      stringsAsFactors = FALSE
    )
  }
}

master <- bind_rows(all_results)

# ---- 9. Display pivots ---------------------------------------------------
strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF")

cat("\n--- MEDIAN RR BY STRATEGY x COMBINATION ---\n")
pivot_median <- master %>%
  select(strategy, combination, RR_median) %>%
  pivot_wider(names_from = combination, values_from = RR_median) %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_median)

cat("\n--- SHORTFALL_70 BY STRATEGY x COMBINATION ---\n")
pivot_sf <- master %>%
  select(strategy, combination, P_shortfall_70) %>%
  pivot_wider(names_from = combination, values_from = P_shortfall_70) %>%
  mutate(strategy = factor(strategy, levels = strategy_order)) %>%
  arrange(strategy)
print(pivot_sf)

# ---- 10. Save ------------------------------------------------------------
saveRDS(master, "Data/robustness_B_regime_thresholds.rds")
cat("\nSaved Data/robustness_B_regime_thresholds.rds\n")
