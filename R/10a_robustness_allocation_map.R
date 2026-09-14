# 10a_robustness_allocation_map.R
#
# Input:  Data/dnb_equity_allocations.rds   (regime predictions are independent
#                                            of the allocation map; only the
#                                            mapping changes here)
#         Data/equity_returns.rds
#         Data/bond_returns_10y.rds
#         Data/salary_path.rds
#         Data/annuity_factor.rds
#         Data/mv_weights.rds
#         Data/ml_models.rds
#         Data/feature_cols.rds
#         Data/features_historical.rds
#         Data/cp2022-...xlsx
#
# Output: Data/robustness_A_alloc_maps.rds  (data.frame: 8 strategies x 4 maps)
#
# Tests the sensitivity of ML strategies' RR distribution to alternative
# regime->equity allocation mappings. Rule-based and MV strategies are
# unaffected. Three alternative maps are compared to baseline.
#
# Baseline (M0):     70 / 55 / 45 / 30  (Normal / Inflationary / Late / Stagflation)
# Aggressive (M1):   80 / 65 / 50 / 35
# Conservative (M2): 60 / 50 / 45 / 40
# Extreme (M3):      80 / 55 / 40 / 20
#
# Note: equity cap from IORP II is 70%, so M1 and M3 (which include 80%)
# violate the cap. They are kept for sensitivity purposes and explicitly
# flagged as "cap-violating" in the output.
#
# Methodology note: this is a PRE-DECLARED sensitivity analysis. Baseline
# results (M0) remain the main thesis numbers.

library(dplyr)
library(glmnet)
library(ranger)
library(readxl)

# ---- 1. Load reusable predictions and inputs -----------------------------
allocations  <- readRDS("Data/dnb_equity_allocations.rds")  # 20000 x 42 x 5
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")
mv_weights   <- readRDS("Data/mv_weights.rds")
models       <- readRDS("Data/ml_models.rds")
feature_cols <- readRDS("Data/feature_cols.rds")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- dim(allocations)[2]
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

# ---- 2. Rebuild ML regime predictions on DNB (not allocations) -----------
# We need the regime predictions per scenario-year so we can apply
# different allocation maps without re-running ML.
#
# To do this, we re-run inference and SAVE the predicted regime class
# (not the equity weight). This is the only way to switch alloc map.

hist <- readRDS("Data/features_historical.rds")
hist_last <- tail(hist, 1)

y1  <- as.matrix(readRDS("Data/yield_1y.rds"))  * 100
y5  <- as.matrix(readRDS("Data/yield_5y.rds"))  * 100
y10 <- as.matrix(readRDS("Data/yield_10y.rds")) * 100
y20 <- as.matrix(readRDS("Data/yield_20y.rds")) * 100
hicp_yoy_dnb <- as.matrix(read_excel(
  "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx",
  sheet = "5_EU_Price_Inflation", col_names = FALSE)) * 100

# Warm-start anchors
warm <- list(
  eq_ret_l1    = hist_last$eq_ret,
  eq_ret_l3    = hist_last$eq_ret_l3,
  eq_ret_l12   = hist_last$eq_ret_l12,
  hicp_yoy_l3  = hist_last$hicp_yoy_l3,
  hicp_yoy_l12 = hist_last$hicp_yoy_l12,
  slope_l3     = hist_last$slope_l3,
  slope_l12    = hist_last$slope_l12,
  y10_l12      = hist_last$y10_l12
)

build_features_year <- function(t_career) {
  t_yc <- t_career + 1
  t_yp <- t_career
  t_rc <- t_career
  t_rp <- max(t_career - 1, 1)
  t_ic <- t_career
  t_ip <- max(t_career - 1, 1)
  
  eq_ret_t   <- eq_ret[, t_rc]
  hicp_yoy_t <- hicp_yoy_dnb[, t_ic]
  y1_t       <- y1[,  t_yc]
  y5_t       <- y5[,  t_yc]
  y10_t      <- y10[, t_yc]
  y20_t      <- y20[, t_yc]
  
  if (t_career == 1) {
    y1_prev  <- y1[, 1]; y10_prev <- y10[, 1]
  } else {
    y1_prev  <- y1[, t_yp]; y10_prev <- y10[, t_yp]
  }
  dy1_t  <- 100 * (y1_t  - y1_prev)
  dy10_t <- 100 * (y10_t - y10_prev)
  slope_10_1_t   <- y10_t - y1_t
  slope_20_5_t   <- y20_t - y5_t
  curvature_t    <- 2 * y5_t - y1_t - y10_t
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
    eq_ret_l1_t   <- eq_ret[, t_rp]
    eq_ret_l3_t   <- eq_ret_l1_t
    eq_ret_l12_t  <- eq_ret[, t_rp]
    hicp_yoy_l3_t  <- hicp_yoy_t
    hicp_yoy_l12_t <- hicp_yoy_dnb[, t_ip]
    y1p <- y1[, t_yp]; y10p <- y10[, t_yp]
    sp  <- y10p - y1p
    slope_l3_t  <- sp
    slope_l12_t <- sp
    y10_l12_t   <- y10p
  }
  slope_x_infl_t <- slope_10_1_t * hicp_yoy_t
  realy_x_eq_t   <- real_yield10_t * eq_ret_l1_t
  
  X <- cbind(
    eq_ret = eq_ret_t, hicp_yoy = hicp_yoy_t,
    y1 = y1_t, y5 = y5_t, y10 = y10_t, y20 = y20_t,
    dy1 = dy1_t, dy10 = dy10_t,
    slope_10_1 = slope_10_1_t, slope_20_5 = slope_20_5_t,
    curvature = curvature_t, real_yield10 = real_yield10_t,
    eq_drawdown_12m = eq_drawdown_12m_t,
    eq_ret_l1 = eq_ret_l1_t, eq_ret_l3 = eq_ret_l3_t, eq_ret_l12 = eq_ret_l12_t,
    hicp_yoy_l3 = hicp_yoy_l3_t, hicp_yoy_l12 = hicp_yoy_l12_t,
    slope_l3 = slope_l3_t, slope_l12 = slope_l12_t, y10_l12 = y10_l12_t,
    slope_x_infl = slope_x_infl_t, realy_x_eq = realy_x_eq_t
  )
  X
}

# ---- 3. Run inference once to get regime CLASS labels per scenario-year --
cat("Running ML inference once to extract regime classes...\n")
ml_classes <- list(
  Lasso = matrix("", nrow = N_SCEN, ncol = N_YEARS),
  Ridge = matrix("", nrow = N_SCEN, ncol = N_YEARS),
  RF    = matrix("", nrow = N_SCEN, ncol = N_YEARS)
)

for (t in 1:N_YEARS) {
  X_t <- build_features_year(t)
  
  pl <- predict(models$lasso, newx = X_t, s = "lambda.min", type = "response")[, , 1]
  ml_classes$Lasso[, t] <- colnames(pl)[apply(pl, 1, which.max)]
  
  pr <- predict(models$ridge, newx = X_t, s = "lambda.min", type = "response")[, , 1]
  ml_classes$Ridge[, t] <- colnames(pr)[apply(pr, 1, which.max)]
  
  pf <- predict(models$rf, data = as.data.frame(X_t))$predictions
  ml_classes$RF[, t] <- colnames(pf)[apply(pf, 1, which.max)]
  
  if (t %% 10 == 0) cat(sprintf("  Year %d / %d\n", t, N_YEARS))
}

# ---- 4. Define alternative allocation maps ------------------------------
maps <- list(
  M0_baseline = c(Normal_expansion = 0.70, Inflationary_exp = 0.55,
                  Late_cycle = 0.45, Stagflation_risk = 0.30),
  M1_aggressive = c(Normal_expansion = 0.80, Inflationary_exp = 0.65,
                    Late_cycle = 0.50, Stagflation_risk = 0.35),
  M2_conservative = c(Normal_expansion = 0.60, Inflationary_exp = 0.50,
                      Late_cycle = 0.45, Stagflation_risk = 0.40),
  M3_extreme = c(Normal_expansion = 0.80, Inflationary_exp = 0.55,
                 Late_cycle = 0.40, Stagflation_risk = 0.20)
)

cat("\n--- ALLOCATION MAPS ---\n")
for (m in names(maps)) {
  cat(m, ":", paste(round(maps[[m]] * 100), collapse = " / "), "%\n")
}

# ---- 5. Wealth and RR helpers --------------------------------------------
compute_wealth_final <- function(w_eq) {
  wealth <- numeric(N_SCEN)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wealth   <- (wealth + contrib) * (1 + port_ret)
  }
  wealth
}

# Rule-based and MV weights stay the same across all maps
w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
w_glide <- matrix(rep(seq(0.80, 0.30, length.out = N_YEARS),
                      each = N_SCEN), N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_weights$plain["equity"]), N_SCEN, N_YEARS)
w_mv_lw <- matrix(as.numeric(mv_weights$lw["equity"]),    N_SCEN, N_YEARS)

# Cache wealth for invariant strategies (computed once)
wf_6040  <- compute_wealth_final(w_6040)
wf_1N    <- compute_wealth_final(w_1N)
wf_glide <- compute_wealth_final(w_glide)
wf_mv_p  <- compute_wealth_final(w_mv_p)
wf_mv_lw <- compute_wealth_final(w_mv_lw)

# ---- 6. Compute RR for each map -----------------------------------------
results <- list()

for (m_name in names(maps)) {
  cat(sprintf("\nProcessing map %s...\n", m_name))
  m <- maps[[m_name]]
  
  # Apply map to ML classes
  alloc_lasso <- matrix(m[ml_classes$Lasso], N_SCEN, N_YEARS)
  alloc_ridge <- matrix(m[ml_classes$Ridge], N_SCEN, N_YEARS)
  alloc_rf    <- matrix(m[ml_classes$RF],    N_SCEN, N_YEARS)
  
  wf_lasso <- compute_wealth_final(alloc_lasso)
  wf_ridge <- compute_wealth_final(alloc_ridge)
  wf_rf    <- compute_wealth_final(alloc_rf)
  
  wealth_all <- cbind(
    `60/40` = wf_6040, `1/N` = wf_1N, Glide = wf_glide,
    MV_plain = wf_mv_p, MV_LW = wf_mv_lw,
    Lasso = wf_lasso, Ridge = wf_ridge, RF = wf_rf
  )
  
  annuity_all <- sweep(wealth_all, 1, af, FUN = "/")
  RR          <- sweep(annuity_all, 1, salary_at_retirement, FUN = "/")
  
  results[[m_name]] <- data.frame(
    strategy        = colnames(RR),
    map             = m_name,
    RR_mean         = round(colMeans(RR), 3),
    RR_median       = round(apply(RR, 2, median), 3),
    RR_sd           = round(apply(RR, 2, sd), 3),
    P_shortfall_70  = round(colMeans(RR < 0.70), 3),
    CVaR_5pct       = round(sapply(colnames(RR), function(s) {
      cutoff <- quantile(RR[, s], 0.05)
      mean(RR[RR[, s] <= cutoff, s])
    }), 3),
    stringsAsFactors = FALSE
  )
}

# ---- 7. Combine into master table ---------------------------------------
master <- bind_rows(results)
cat("\n--- MASTER ROBUSTNESS-A TABLE ---\n\n")
print(master, row.names = FALSE)

# Compare: how does each strategy's median RR change across maps?
cat("\n--- MEDIAN RR BY STRATEGY x MAP ---\n")
pivot <- master %>%
  select(strategy, map, RR_median) %>%
  pivot_wider(names_from = map, values_from = RR_median)
print(pivot)

cat("\n--- SHORTFALL_70 BY STRATEGY x MAP ---\n")
pivot_sf <- master %>%
  select(strategy, map, P_shortfall_70) %>%
  pivot_wider(names_from = map, values_from = P_shortfall_70)
print(pivot_sf)

# ---- 8. Save -------------------------------------------------------------
saveRDS(master, "Data/robustness_A_alloc_maps.rds")
cat("\nSaved Data/robustness_A_alloc_maps.rds\n")

# Reminder: maps M1 and M3 contain a 0.80 weight, which violates the
# IORP II 70% equity cap. Document this explicitly in the thesis.
cat("\nNOTE: M1 and M3 exceed IORP II 70% equity cap. Reported for sensitivity.\n")
