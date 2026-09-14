# EXPORT_RESULTS.R 
# ============================================================

cat("################################################################\n")
cat("# SECTION 5 RESULTS — DATA EXPORT\n")
cat("################################################################\n\n")

# ---- 1. MASTER METRICS TABLE -------------------------------------
cat("============================================================\n")
cat("1. MASTER METRICS TABLE (all strategies)\n")
cat("============================================================\n")
metrics <- readRDS("Data/metrics_full.rds")

#  All columns, full format
cat("\n--- FULL TABLE (all columns) ---\n")
print(metrics, row.names = FALSE, width = 200)

cat("\n--- PERFORMANCE COLUMNS ---\n")
perf_cols <- c("strategy", "RR_mean", "RR_median", "RR_sd", 
               "RR_p10", "RR_p25", "RR_p75", "RR_p90", "RR_iqr", "RR_sharpe_like")
perf_cols <- intersect(perf_cols, colnames(metrics))
print(metrics[, perf_cols], row.names = FALSE)

cat("\n--- RISK COLUMNS ---\n")
risk_cols <- c("strategy", "P_shortfall_50pct", "P_shortfall_60pct", 
               "P_shortfall_70pct", "P_shortfall_80pct",
               "CVaR_5pct", "CVaR_10pct", "MDD_mean", "MDD_median", "MDD_p10",
               "turnover_mean", "turnover_median")
risk_cols <- intersect(risk_cols, colnames(metrics))
print(metrics[, risk_cols], row.names = FALSE)

# ---- 2. RANKINGS -------------------------------------------------
cat("\n\n============================================================\n")
cat("2. RANKINGS BY KEY METRIC\n")
cat("============================================================\n")

cat("\n--- By RR_median (desc) ---\n")
o <- order(-metrics$RR_median)
print(data.frame(rank = 1:nrow(metrics),
                 strategy = metrics$strategy[o],
                 RR_median = metrics$RR_median[o]), row.names = FALSE)

cat("\n--- By P_shortfall_70pct (asc, lower=better) ---\n")
o <- order(metrics$P_shortfall_70pct)
print(data.frame(rank = 1:nrow(metrics),
                 strategy = metrics$strategy[o],
                 shortfall_70 = metrics$P_shortfall_70pct[o]), row.names = FALSE)

cat("\n--- By CVaR_5pct (desc, higher=better) ---\n")
o <- order(-metrics$CVaR_5pct)
print(data.frame(rank = 1:nrow(metrics),
                 strategy = metrics$strategy[o],
                 CVaR_5 = metrics$CVaR_5pct[o]), row.names = FALSE)

cat("\n--- By MDD_mean (desc, less negative=better) ---\n")
o <- order(-metrics$MDD_mean)
print(data.frame(rank = 1:nrow(metrics),
                 strategy = metrics$strategy[o],
                 MDD_mean = metrics$MDD_mean[o]), row.names = FALSE)

cat("\n--- By turnover_mean (asc, lower=better) ---\n")
o <- order(metrics$turnover_mean)
print(data.frame(rank = 1:nrow(metrics),
                 strategy = metrics$strategy[o],
                 turnover = metrics$turnover_mean[o]), row.names = FALSE)

# ---- 3. RR DISTRIBUTION DETAILS ----------------------------------
cat("\n\n============================================================\n")
cat("3. RR DISTRIBUTION (extra percentiles for figures)\n")
cat("============================================================\n")
RR <- tryCatch(readRDS("Data/replacement_ratios_with_hybrid_naive.rds"),
               error = function(e) NULL)
if (!is.null(RR)) {
  cat("RR matrix dim:", dim(RR), "\n")
  cat("Columns:", paste(colnames(RR), collapse=", "), "\n\n")
  cat("--- Percentiles per strategy (p1, p5, p10, p25, p50, p75, p90, p95, p99) ---\n")
  probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  pct_table <- t(apply(RR, 2, quantile, probs = probs))
  colnames(pct_table) <- paste0("p", probs*100)
  print(round(pct_table, 4))
  
  cat("\n--- Mean, SD, Skewness, Min, Max ---\n")
  skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
  summ <- data.frame(
    strategy = colnames(RR),
    mean = round(apply(RR, 2, mean), 4),
    sd = round(apply(RR, 2, sd), 4),
    skewness = round(apply(RR, 2, skew), 4),
    min = round(apply(RR, 2, min), 4),
    max = round(apply(RR, 2, max), 4)
  )
  print(summ, row.names = FALSE)
} else {
  cat("File not found, trying replacement_ratios_all.rds\n")
  RR <- readRDS("Data/replacement_ratios_all.rds")
  cat("RR matrix dim:", dim(RR), "\n")
  cat("Columns:", paste(colnames(RR), collapse=", "), "\n")
}

# ---- 4. ROBUSTNESS CHECKS ----------------------------------------
cat("\n\n============================================================\n")
cat("4. ROBUSTNESS CHECKS\n")
cat("============================================================\n")

rob_files <- list(
  "10a_allocation_map"   = "Data/robustness_A_alloc_maps.rds",
  "10b_regime_threshold" = "Data/robustness_B_regime_thresholds.rds",
  "10c_drawdown"         = "Data/robustness_C_drawdown_threshold.rds",
  "10d_bond_duration"    = "Data/robustness_D_bond_duration.rds",
  "10e_mv_window"        = "Data/robustness_E_mv_window.rds",
  "10h_gamma"            = "Data/robustness_H_gamma.rds",
  "10i_contribution"     = "Data/robustness_I_contribution.rds",
  "10k_equity_cap"       = "Data/robustness_K_equity_cap.rds",
  "10L_target"           = "Data/robustness_L_target.rds",
  "10q_wage_drift"       = "Data/robustness_Q_wage_drift.rds"
)

for (nm in names(rob_files)) {
  f <- rob_files[[nm]]
  cat("\n------------------------------------------------------------\n")
  cat("CHECK:", nm, " (", f, ")\n")
  cat("------------------------------------------------------------\n")
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) {
    cat("  [FILE NOT FOUND]\n")
    next
  }
  cat("  Class:", class(obj), "\n")
  if (is.data.frame(obj)) {
    cat("  Dim:", dim(obj), "\n")
    print(obj, row.names = FALSE, width = 200)
  } else if (is.list(obj)) {
    cat("  List names:", paste(names(obj), collapse=", "), "\n")
    for (el in names(obj)) {
      cat("\n  >> Element:", el, "\n")
      x <- obj[[el]]
      if (is.data.frame(x)) {
        print(x, row.names = FALSE, width = 200)
      } else if (is.matrix(x)) {
        print(round(x, 4))
      } else {
        print(x)
      }
    }
  } else if (is.matrix(obj)) {
    print(round(obj, 4))
  } else {
    print(obj)
  }
}

cat("\n\n################################################################\n")
cat("# EXPORT COMPLETE\n")
cat("################################################################\n")


N_YEARS <- 42

eq_ret_full <- readRDS("Data/equity_returns.rds")        # 20000 x 100
bd_ret_full <- readRDS("Data/bond_returns_10y.rds")      # 20000 x 100
allocations <- readRDS("Data/dnb_equity_allocations.rds") # 20000 x 42 x 5

# Cut to first 42 years (accumulation horizon)
eq_ret <- eq_ret_full[, 1:N_YEARS]
bd_ret <- bd_ret_full[, 1:N_YEARS]
N_SCEN <- nrow(eq_ret)

cat("Dimension: eq_ret", dim(eq_ret)[1], "x", dim(eq_ret)[2], "\n")
cat("Strategy in allocations:", dimnames(allocations)[[3]], "\n\n")


w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
glide_v <- seq(0.80, 0.30, length.out = N_YEARS)
w_glide <- matrix(rep(glide_v, each = N_SCEN), N_SCEN, N_YEARS)
w_mvp   <- allocations[, , "MV_plain"]
w_mvlw  <- allocations[, , "MV_LW"]
w_lasso <- allocations[, , "Lasso"]
w_ridge <- allocations[, , "Ridge"]
w_rf    <- allocations[, , "RF"]

# Hybrid: average 60/40, MV_plain, RF 
w_hyb <- (w_6040 + w_mvp + w_rf) / 3

strat <- list("60/40"=w_6040, "1/N"=w_1N, "Glide"=w_glide,
              "MV_plain"=w_mvp, "MV_LW"=w_mvlw,
              "Lasso"=w_lasso, "Ridge"=w_ridge, "RF"=w_rf,
              "HybridNaive"=w_hyb)

# --- Time-weighted return ---
twr_annual <- function(w_eq) {
  port_ret <- w_eq * eq_ret + (1 - w_eq) * bd_ret
  growth   <- apply(1 + port_ret, 1, prod)
  growth^(1/N_YEARS) - 1
}

results <- data.frame(strategy=names(strat),
                      median_ann=NA, mean_ann=NA, p10=NA, p90=NA)
for (i in seq_along(strat)) {
  ann <- twr_annual(strat[[i]])
  results$median_ann[i] <- round(median(ann)*100, 2)
  results$mean_ann[i]   <- round(mean(ann)*100, 2)
  results$p10[i]        <- round(quantile(ann,0.10)*100, 2)
  results$p90[i]        <- round(quantile(ann,0.90)*100, 2)
}
results <- results[order(-results$median_ann), ]
cat("Annualized time-weighted return (%, nominal):\n\n")
print(results, row.names=FALSE)

cat("\n--- Mean equity weight (context) ---\n")
for (i in seq_along(strat))
  cat(sprintf("  %-12s : %.3f\n", names(strat)[i], mean(strat[[i]])))

# Cosa contiene la lista 'rr'
cat("=== elementi di rr ===\n"); print(names(rr))
