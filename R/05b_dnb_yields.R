# 05b_dnb_yields.R
#
# Input:        Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx
# Output:       Data/yield_1y.rds   (20000 x 101)
#               Data/yield_5y.rds   (20000 x 101)
#               Data/yield_10y.rds  (20000 x 101)  -- already exists, re-verified
#               Data/yield_20y.rds  (20000 x 101)  -- already exists, re-verified
# Dependencies: readxl
#
# Constructs DNB-derived yields at maturities 1y, 5y, 10y, 20y for all
# 20,000 scenarios and 101 time points (t=0 to t=100) using the
# Commissie Parameters 2022 nominal yield formula:
#
#   y(t, tau) = exp(-bracket / tau) - 1
#
# where:
#   bracket = Phi_t(tau) + x1_t * Psi_1(tau) + x2_t * Psi_2(tau) + x3_t * Psi_3(tau)
#
# Phi:   maturities x time matrix from sheet "7_Nominal_interest_rate_par_Phi"
# Psi:   maturities x state_vars matrix from sheet "8_Nominal_interest_rate_par_Psi"
# x1,2,3: state variables from sheets "1_..", "2_..", "3_..".

library(readxl)

file_path <- "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx"

# ---- 1. Load all required inputs -----------------------------------------
cat("Loading DNB data sheets...\n")
state_var_1 <- read_excel(file_path, sheet = "1_Statevariable_1", col_names = FALSE)
state_var_2 <- read_excel(file_path, sheet = "2_Statevariable_2", col_names = FALSE)
state_var_3 <- read_excel(file_path, sheet = "3_Statevariable_3", col_names = FALSE)
phi_params  <- read_excel(file_path, sheet = "7_Nominal_interest_rate_par_Phi", col_names = FALSE)
psi_params  <- read_excel(file_path, sheet = "8_Nominal_interest_rate_par_Psi", col_names = FALSE)

# Convert to matrices
phi_mat <- as.matrix(phi_params)   # rows = maturities (1..100), cols = time (0..100)
psi_mat <- as.matrix(psi_params)   # rows = maturities (1..100), cols = state vars (1..3)
x1_mat  <- as.matrix(state_var_1)  # rows = scenarios (20000), cols = time (0..100)
x2_mat  <- as.matrix(state_var_2)
x3_mat  <- as.matrix(state_var_3)

cat("Phi  dim:", dim(phi_mat), "\n")
cat("Psi  dim:", dim(psi_mat), "\n")
cat("x1   dim:", dim(x1_mat),  "\n\n")

# ---- 2. Vectorised yield-curve computation -------------------------------
# For a SINGLE target maturity tau, the formula computes one yield value
# per (scenario, time). We can vectorise over all (scenario, time) pairs:
#
#   For each time t and maturity tau (fixed):
#     bracket[i, t] = Phi[tau, t] + x1[i, t]*Psi[tau, 1]
#                   + x2[i, t]*Psi[tau, 2] + x3[i, t]*Psi[tau, 3]
#     yield[i, t]   = exp(-bracket[i, t] / tau) - 1

compute_yields_all_scenarios <- function(target_maturity) {
  n_scenarios <- nrow(x1_mat)
  n_time      <- ncol(x1_mat)  # 101
  tau         <- target_maturity
  
  # Phi at this maturity for each time t: vector of length n_time
  phi_tau <- phi_mat[tau, ]            # 1 x n_time
  
  # Psi at this maturity: 3 values
  psi_tau1 <- psi_mat[tau, 1]
  psi_tau2 <- psi_mat[tau, 2]
  psi_tau3 <- psi_mat[tau, 3]
  
  # Broadcast phi_tau over scenarios (replicate to n_scenarios rows)
  phi_broadcast <- matrix(phi_tau, nrow = n_scenarios, ncol = n_time, byrow = TRUE)
  
  # Bracket: vectorised over all (scenario, time)
  bracket <- phi_broadcast +
    x1_mat * psi_tau1 +
    x2_mat * psi_tau2 +
    x3_mat * psi_tau3
  
  # Yield
  yields <- exp(-bracket / tau) - 1
  yields
}

# ---- 3. Compute yields for 1y, 5y, 10y, 20y ------------------------------
maturities <- c(1, 5, 10, 20)

for (tau in maturities) {
  cat(sprintf("Computing yield curve at maturity %dy...\n", tau))
  start_time <- Sys.time()
  yields_mat <- compute_yields_all_scenarios(tau)
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  cat(sprintf("  Done in %.1fs. Dim: %d x %d\n",
              elapsed, nrow(yields_mat), ncol(yields_mat)))
  
  # Sanity checks
  cat(sprintf("  Range:      [%.4f, %.4f]\n", min(yields_mat), max(yields_mat)))
  cat(sprintf("  Mean:        %.4f\n", mean(yields_mat)))
  cat(sprintf("  t=0 value:   %.6f (should be identical across all scenarios)\n",
              yields_mat[1, 1]))
  cat(sprintf("  t=0 unique:  %d (should be 1)\n\n",
              length(unique(yields_mat[, 1]))))
  
  # Save
  out_path <- sprintf("Data/yield_%dy.rds", tau)
  saveRDS(yields_mat, out_path)
  cat(sprintf("  Saved %s\n\n", out_path))
}

# ---- 4. Cross-maturity consistency check ---------------------------------
cat("--- CROSS-MATURITY CONSISTENCY (scenario 1, first 5 time points) ---\n")
y1  <- readRDS("Data/yield_1y.rds")
y5  <- readRDS("Data/yield_5y.rds")
y10 <- readRDS("Data/yield_10y.rds")
y20 <- readRDS("Data/yield_20y.rds")

check <- data.frame(
  t   = 0:4,
  y1  = y1 [1, 1:5],
  y5  = y5 [1, 1:5],
  y10 = y10[1, 1:5],
  y20 = y20[1, 1:5]
)
print(round(check, 4))

cat("\n--- SLOPE 10y-1y (scenario 1, first 5 time points) ---\n")
print(round(y10[1, 1:5] - y1[1, 1:5], 4))

cat("\nAll yields ready for use in 06_apply_to_dnb.R\n")
