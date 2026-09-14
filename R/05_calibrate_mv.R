# 05_calibrate_mv.R
#
# Input:        Data/features_historical.rds
# Output:       Data/mv_weights.rds       (production MV weights)
#               Data/mv_calibration.rds   (calibration moments + diagnostics)
# Dependencies: dplyr, lubridate
#
# Calibrates Mean-Variance optimal portfolio weights on full historical sample
# (2001-2024, 288 monthly observations). Two variants:
#   - MV plain:        sample covariance estimator
#   - MV Ledoit-Wolf:  shrinkage covariance estimator toward identity target
#
# Static calibration: a single (mu, Sigma) is estimated, yielding a single
# optimal weight vector w* applied to all DNB scenarios. This mirrors the
# pension fund advisor's strategic allocation framework, where MV is used
# as a calibration tool, not a tactical signal.
#
# Setup:
#   Assets: equity (MSCI World NTR EUR), bond (rolling 10y Bund total return)
#   Risk aversion gamma = 5 (Campbell & Viceira 2002)
#   Constraints: w >= 0, sum(w) = 1, w_equity <= 0.70 (IORP II Article 19)

library(dplyr)
library(lubridate)

# ---- 1. Load and build asset returns -------------------------------------
df <- readRDS("Data/features_historical.rds") %>%
  arrange(date) %>%
  select(date, eq_price, y10) %>%
  mutate(
    # Equity log return (already in features set, recomputed here for clarity)
    eq_ret = log(eq_price / lag(eq_price, 1)),
    # Bond return approximation: total return of a constant-maturity 10y bond
    # held over one month, given month-on-month yield change.
    # Formula: r_bond_t ~ y10_{t-1}/12  - D * (y10_t - y10_{t-1}) / 100
    # where D = modified duration ~ 9 (for a 10y bond at typical yields).
    # Yields are in % so the change must be in decimal form.
    D = 9,
    y10_lag = lag(y10, 1),
    bond_ret = (y10_lag / 100) / 12 - D * (y10 - y10_lag) / 100
  ) %>%
  select(-D, -y10_lag) %>%
  filter(!is.na(eq_ret) & !is.na(bond_ret))

cat("--- ASSET RETURN SUMMARY (monthly) ---\n")
ret_mat <- as.matrix(df[, c("eq_ret", "bond_ret")])
print(summary(ret_mat))

cat("\nN observations:", nrow(ret_mat), "\n")
cat("Date range:    ", as.character(min(df$date)),
    "to", as.character(max(df$date)), "\n")

# ---- 2. Annualised moments -----------------------------------------------
# Convert monthly moments to annual for MV (since DNB is annual and gamma
# is calibrated at annual horizon).
mu_m  <- colMeans(ret_mat)
Sig_m <- cov(ret_mat)

mu_a  <- mu_m * 12               # annualised mean
Sig_a <- Sig_m * 12              # annualised covariance (i.i.d. assumption)

cat("\n--- ANNUALISED MOMENTS ---\n")
cat("Mean returns (annual):\n")
print(round(mu_a, 4))
cat("\nCovariance matrix (annual):\n")
print(round(Sig_a, 6))
cat("\nVolatilities (annual):\n")
print(round(sqrt(diag(Sig_a)), 4))
cat("\nCorrelation:\n")
print(round(cov2cor(Sig_a), 4))

# ---- 3. Ledoit-Wolf shrinkage covariance ---------------------------------
# Shrinks sample covariance toward the identity target:
#   Sig_lw = (1 - delta) * Sig_sample + delta * F
# where F is the structured target (scaled identity) and delta is the
# optimal shrinkage intensity (Ledoit & Wolf 2003 closed-form).
#
# For 2 assets, LW shrinkage is mild (sample cov is already low-dimensional).
# This is still useful to demonstrate the methodology and provides a
# robustness comparison with the plain estimator.

ledoit_wolf_shrinkage <- function(R) {
  # R: T x N return matrix
  T <- nrow(R); N <- ncol(R)
  R_centered <- scale(R, center = TRUE, scale = FALSE)
  S <- cov(R)  # sample covariance
  # Target: scaled identity (mean of diagonal)
  mu_diag <- mean(diag(S))
  F <- diag(mu_diag, N)
  # Shrinkage intensity (simplified closed-form)
  # pi: sum of asymptotic variances of S entries
  pi_mat <- (t(R_centered^2) %*% (R_centered^2)) / T - S^2
  pi_hat <- sum(pi_mat)
  # rho: asymptotic covariance term (zero for identity target with same diag)
  rho_hat <- sum(diag(pi_mat))
  # gamma: distance between sample and target
  gamma_hat <- sum((S - F)^2)
  # delta optimal
  kappa <- (pi_hat - rho_hat) / gamma_hat
  delta <- max(0, min(1, kappa / T))
  S_lw <- delta * F + (1 - delta) * S
  list(Sigma = S_lw, delta = delta, S_sample = S, F = F)
}

lw_result <- ledoit_wolf_shrinkage(ret_mat)
Sig_lw_m  <- lw_result$Sigma
Sig_lw_a  <- Sig_lw_m * 12  # annualise
delta_lw  <- lw_result$delta

cat("\n--- LEDOIT-WOLF SHRINKAGE ---\n")
cat(sprintf("Shrinkage intensity delta: %.4f\n", delta_lw))
cat("LW-shrunk covariance (annual):\n")
print(round(Sig_lw_a, 6))

# ---- 4. MV optimisation with constraints ---------------------------------
# Solve: max mu' w - (gamma/2) * w' Sigma w
# s.t.   sum(w) = 1, w >= 0, w_equity <= 0.70
#
# With 2 assets and constraints, the problem is one-dimensional in w_equity.
# We grid-search w_equity in [0, 0.70] and pick the maximiser. This is
# precise enough for 2 assets and avoids dependency on quadprog/CVXR.

GAMMA <- 5
EQ_CAP <- 0.70

mv_optimal_weights <- function(mu, Sigma, gamma = 5, eq_cap = 0.70) {
  # mu, Sigma: 2-vector and 2x2 matrix (assets ordered: equity, bond)
  w_eq_grid <- seq(0, eq_cap, by = 0.001)
  utility <- sapply(w_eq_grid, function(we) {
    w <- c(we, 1 - we)
    sum(w * mu) - (gamma / 2) * t(w) %*% Sigma %*% w
  })
  w_eq_opt <- w_eq_grid[which.max(utility)]
  c(equity = w_eq_opt, bond = 1 - w_eq_opt)
}

w_plain <- mv_optimal_weights(mu_a, Sig_a,    gamma = GAMMA, eq_cap = EQ_CAP)
w_lw    <- mv_optimal_weights(mu_a, Sig_lw_a, gamma = GAMMA, eq_cap = EQ_CAP)

cat("\n--- OPTIMAL MV WEIGHTS ---\n")
cat("MV plain:        equity =", round(w_plain["equity"], 3),
    ", bond =", round(w_plain["bond"], 3), "\n")
cat("MV Ledoit-Wolf:  equity =", round(w_lw["equity"], 3),
    ", bond =", round(w_lw["bond"], 3), "\n")

# Flag if equity cap binds (constraint is active)
if (w_plain["equity"] >= EQ_CAP - 0.001) {
  cat("\nNote: equity cap binding for MV plain (interior solution would prefer higher equity)\n")
}
if (w_lw["equity"] >= EQ_CAP - 0.001) {
  cat("Note: equity cap binding for MV Ledoit-Wolf\n")
}

# ---- 5. Sanity check: implied portfolio statistics -----------------------
ptf_stats <- function(w, mu, Sigma) {
  list(
    expected_return = sum(w * mu),
    volatility      = sqrt(t(w) %*% Sigma %*% w),
    sharpe          = sum(w * mu) / sqrt(t(w) %*% Sigma %*% w)
  )
}

cat("\n--- IMPLIED PORTFOLIO STATISTICS (annualised) ---\n")
stats_plain <- ptf_stats(w_plain, mu_a, Sig_a)
stats_lw    <- ptf_stats(w_lw, mu_a, Sig_lw_a)

cat(sprintf("MV plain:        E[R] = %.3f, vol = %.3f, SR = %.3f\n",
            stats_plain$expected_return, stats_plain$volatility, stats_plain$sharpe))
cat(sprintf("MV Ledoit-Wolf:  E[R] = %.3f, vol = %.3f, SR = %.3f\n",
            stats_lw$expected_return, stats_lw$volatility, stats_lw$sharpe))

# ---- 6. Robustness: alternative calibration windows ----------------------
# Pre-declared in decisions_log.md: window in {60, 120, 240} months.
# Default uses full sample (288 months).
cat("\n--- ROBUSTNESS: ALTERNATIVE CALIBRATION WINDOWS ---\n")

windows <- c(60, 120, 240, nrow(ret_mat))
names(windows) <- c("60m", "120m", "240m", "full_sample")

robustness <- data.frame(
  window      = names(windows),
  n_obs       = numeric(length(windows)),
  w_eq_plain  = numeric(length(windows)),
  w_eq_lw     = numeric(length(windows)),
  delta_lw    = numeric(length(windows))
)

for (i in seq_along(windows)) {
  w_len <- windows[i]
  R_win <- ret_mat[(nrow(ret_mat) - w_len + 1):nrow(ret_mat), ]
  mu_win <- colMeans(R_win) * 12
  Sig_win <- cov(R_win) * 12
  lw_win <- ledoit_wolf_shrinkage(R_win)
  Sig_lw_win <- lw_win$Sigma * 12
  
  w_p <- mv_optimal_weights(mu_win, Sig_win, GAMMA, EQ_CAP)
  w_l <- mv_optimal_weights(mu_win, Sig_lw_win, GAMMA, EQ_CAP)
  
  robustness$n_obs[i]      <- w_len
  robustness$w_eq_plain[i] <- round(w_p["equity"], 3)
  robustness$w_eq_lw[i]    <- round(w_l["equity"], 3)
  robustness$delta_lw[i]   <- round(lw_win$delta, 4)
}

print(robustness)

# ---- 7. Save -------------------------------------------------------------
mv_calibration <- list(
  mu_annual    = mu_a,
  Sigma_plain  = Sig_a,
  Sigma_lw     = Sig_lw_a,
  delta_lw     = delta_lw,
  gamma        = GAMMA,
  equity_cap   = EQ_CAP,
  n_obs        = nrow(ret_mat),
  date_first   = min(df$date),
  date_last    = max(df$date),
  robustness   = robustness
)

mv_weights <- list(
  plain = w_plain,
  lw    = w_lw
)

saveRDS(mv_weights, "Data/mv_weights.rds")
saveRDS(mv_calibration, "Data/mv_calibration.rds")

cat("\nSaved Data/mv_weights.rds\n")
cat("Saved Data/mv_calibration.rds\n")

library(readxl)

dnb_path <- "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx"

# 1. List all sheets
sheets <- excel_sheets(dnb_path)
cat("--- DNB FILE SHEETS ---\n")
print(sheets)

# 2. For each sheet, defensive audit (skip empty)
for (s in sheets) {
  cat("\n\n--- SHEET:", s, "---\n")
  preview <- tryCatch(
    read_excel(dnb_path, sheet = s, n_max = 5),
    error = function(e) { cat("ERROR reading sheet:", e$message, "\n"); NULL }
  )
  if (is.null(preview)) next
  
  if (ncol(preview) == 0) {
    cat("Sheet is empty (0 columns).\n")
    next
  }
  
  cat("Dim (first 5 rows shown):", dim(preview), "\n")
  cat("Column names (first 15):\n")
  print(head(colnames(preview), 15))
  cat("Preview (first ", min(6, ncol(preview)), " columns):\n", sep = "")
  print(preview[, 1:min(6, ncol(preview))])
}
