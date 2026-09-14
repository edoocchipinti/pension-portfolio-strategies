# 13_unsupervised_regimes.R  (V2 — robust EM with standardization + multi-seed retry)
#
# Unsupervised regime detection via 4-state Gaussian HMM.
#
# V2 CHANGES vs V1:
#   - Features are standardized (mean 0, sd 1) before HMM fitting. EM is
#     sensitive to feature scale; raw scales (slope ~1, infl ~2,
#     drawdown ~-0.06) cause convergence to invalid likelihoods.
#   - State means/sds back-transformed to raw scale after fit (for
#     interpretability and for DNB Viterbi which uses raw DNB features).
#   - Multi-seed retry: if EM fails, try 5 different seeds before giving up.
#
# Setup (unchanged):
#   - K = 4 states (matches V3)
#   - 4-dim Gaussian emissions
#   - BIC diagnostic k = 2..6
#   - Hard mapping (Viterbi) + Soft mapping (forward-backward posteriors)
#   - DNB inference: Viterbi+FB on each scenario's 42-year sequence
#
# Dependencies: depmixS4, matrixStats, readxl, tidyr, dplyr

library(depmixS4)
library(matrixStats)
library(readxl)
library(tidyr)
library(dplyr)

# Force-bind dplyr verbs to prevent namespace masking by MASS/nnet (loaded
# by depmixS4)
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise
bind_rows <- dplyr::bind_rows

set.seed(42)

# ---- 1. Load inputs ------------------------------------------------------
hist         <- readRDS("Data/features_historical.rds")
allocations  <- readRDS("Data/dnb_equity_allocations.rds")
eq_ret       <- readRDS("Data/equity_returns.rds")
bd_ret       <- readRDS("Data/bond_returns_10y.rds")
salary       <- readRDS("Data/salary_path.rds")
af           <- readRDS("Data/annuity_factor.rds")

y1  <- as.matrix(readRDS("Data/yield_1y.rds"))  * 100
y5  <- as.matrix(readRDS("Data/yield_5y.rds"))  * 100
y10 <- as.matrix(readRDS("Data/yield_10y.rds")) * 100
y20 <- as.matrix(readRDS("Data/yield_20y.rds")) * 100
hicp_yoy_dnb <- as.matrix(read_excel(
  "Data/cp2022-p-scenarioset-20k-2026q1 (ENG).xlsx",
  sheet = "5_EU_Price_Inflation", col_names = FALSE)) * 100

metrics_old <- readRDS("Data/metrics_full.rds")
wp_old      <- readRDS("Data/wealth_paths.rds")
mdd_old     <- readRDS("Data/max_drawdowns.rds")
to_old      <- readRDS("Data/turnovers.rds")
RR_old      <- readRDS("Data/replacement_ratios_with_hybrid_naive.rds")

# Idempotency: strip previous HMM rows if present
for (s in c("HMM_hard", "HMM_soft")) {
  if (s %in% metrics_old$strategy) {
    cat(sprintf("Removing previous %s rows for re-run\n", s))
    metrics_old <- metrics_old %>% filter(strategy != s)
    wp_old[[s]] <- NULL
    mdd_old <- mdd_old[, colnames(mdd_old) != s, drop = FALSE]
    to_old  <- to_old [, colnames(to_old)  != s, drop = FALSE]
    RR_old  <- RR_old [, colnames(RR_old)  != s, drop = FALSE]
  }
}
cat("\n")

N_SCEN  <- nrow(eq_ret)
N_YEARS <- 42
CONTRIB_RATE <- 0.14
salary_at_retirement <- salary[, N_YEARS + 1]

# ---- 2. Build clustering feature matrix and STANDARDIZE -----------------
clust_features <- c("slope_10_1", "hicp_yoy", "eq_drawdown_12m", "real_yield10")

hist <- hist[order(hist$date), ]
hist_c <- hist[, c("date", clust_features), drop = FALSE]
n_hist <- nrow(hist_c)

cat("Clustering features (historical, RAW scale):\n")
print(summary(hist_c[, clust_features]))

# Compute scaling parameters (store for DNB inference back-transform)
X_raw <- as.matrix(hist_c[, clust_features])
feat_mean <- colMeans(X_raw)
feat_sd   <- apply(X_raw, 2, sd)

cat("\nScaling parameters (mean / sd per feature):\n")
print(round(rbind(mean = feat_mean, sd = feat_sd), 3))

# Standardize
X_std <- scale(X_raw, center = feat_mean, scale = feat_sd)
X_clust <- as.data.frame(X_std)
colnames(X_clust) <- clust_features

cat("\nClustering features (STANDARDIZED, used for EM):\n")
print(summary(X_clust))
cat(sprintf("\nN observations: %d\n\n", n_hist))

# ---- 3. Robust HMM fit with multi-seed retry ----------------------------
# Try multiple seeds; if EM fails on one, try next. Return first
# successful fit.

fit_hmm_robust <- function(k, data, features, seeds = c(42, 7, 123, 2024, 999, 555),
                           maxit = 500, tol = 1e-5) {
  rsp_list <- lapply(features, function(f) {
    as.formula(paste0(f, " ~ 1"))
  })
  fam_list <- lapply(features, function(f) gaussian())
  
  best_fit <- NULL
  best_ll  <- -Inf
  
  for (sd_val in seeds) {
    set.seed(sd_val)
    mod <- depmix(response = rsp_list,
                  data = data,
                  nstates = k,
                  family = fam_list)
    
    fit_obj <- tryCatch(
      depmixS4::fit(mod, verbose = FALSE,
                    emcontrol = em.control(maxit = maxit, tol = tol)),
      error = function(e) {
        cat(sprintf("    Seed %d failed: %s\n", sd_val,
                    substr(e$message, 1, 80)))
        NULL
      }
    )
    
    if (!is.null(fit_obj)) {
      ll <- as.numeric(logLik(fit_obj))
      if (is.finite(ll) && ll > best_ll) {
        best_fit <- fit_obj
        best_ll  <- ll
      }
    }
  }
  
  if (is.null(best_fit)) {
    cat(sprintf("    All seeds failed for k = %d\n", k))
  } else {
    cat(sprintf("    Best logLik across seeds: %.2f\n", best_ll))
  }
  best_fit
}

# ---- 4. BIC diagnostic for k = 2..6 -------------------------------------
cat("Running BIC diagnostic for k = 2..6 (with multi-seed retry)...\n")
bic_results <- data.frame(k = 2:6, BIC = NA_real_, logLik = NA_real_,
                          converged = FALSE)
for (i in seq_along(bic_results$k)) {
  k_val <- bic_results$k[i]
  cat(sprintf("  Fitting k = %d...\n", k_val))
  fit_k <- fit_hmm_robust(k_val, X_clust, clust_features)
  if (!is.null(fit_k)) {
    bic_results$BIC[i]    <- BIC(fit_k)
    bic_results$logLik[i] <- as.numeric(logLik(fit_k))
    bic_results$converged[i] <- TRUE
  }
}

cat("\n--- BIC DIAGNOSTIC ---\n")
print(bic_results)
if (all(is.na(bic_results$BIC))) {
  stop("All HMM fits failed across all k. Investigate data or seed issues.")
}
best_k_bic <- bic_results$k[which.min(bic_results$BIC)]
cat(sprintf("\nBest k by BIC: %d  (we use k=4 for V3 comparability)\n\n",
            best_k_bic))

# ---- 5. Main fit with k = 4 ---------------------------------------------
K <- 4
cat(sprintf("Fitting main HMM with k = %d on %d obs (standardized)...\n",
            K, n_hist))
fit_main <- fit_hmm_robust(K, X_clust, clust_features)
if (is.null(fit_main)) stop("Main HMM fit failed across all seeds.")

cat("Convergence message: ", fit_main@message, "\n")
cat(sprintf("Final log-likelihood: %.2f\n", as.numeric(logLik(fit_main))))
cat(sprintf("BIC: %.2f\n\n", BIC(fit_main)))

# ---- 6. Extract state-specific Gaussian params (in STANDARDIZED scale) -
state_means_std <- matrix(NA_real_, nrow = K, ncol = length(clust_features))
colnames(state_means_std) <- clust_features
rownames(state_means_std) <- paste0("State_", 1:K)
state_sds_std <- state_means_std

for (s in 1:K) {
  for (j in seq_along(clust_features)) {
    state_means_std[s, j] <- fit_main@response[[s]][[j]]@parameters$coefficients
    state_sds_std[s, j]   <- fit_main@response[[s]][[j]]@parameters$sd
  }
}

# Back-transform to RAW scale (for interpretability and for DNB Viterbi)
state_means <- sweep(state_means_std, 2, feat_sd, FUN = "*")
state_means <- sweep(state_means,     2, feat_mean, FUN = "+")
state_sds   <- sweep(state_sds_std,   2, feat_sd, FUN = "*")
colnames(state_means) <- clust_features
colnames(state_sds)   <- clust_features
rownames(state_means) <- paste0("State_", 1:K)
rownames(state_sds)   <- paste0("State_", 1:K)

cat("--- STATE MEAN PARAMETERS (raw scale) ---\n")
print(round(state_means, 3))
cat("\n--- STATE SD PARAMETERS (raw scale) ---\n")
print(round(state_sds, 3))

# Transition matrix
trans_mat <- matrix(NA_real_, K, K)
for (i in 1:K) {
  trans_mat[i, ] <- fit_main@transition[[i]]@parameters$coefficients
}
cat("\n--- TRANSITION MATRIX ---\n")
print(round(trans_mat, 3))
cat("(Row = from state, Column = to state)\n")
cat(sprintf("Mean diagonal probability (persistence): %.3f\n\n",
            mean(diag(trans_mat))))

# ---- 7. Map states to economic regime labels ---------------------------
state_score <- data.frame(
  state = 1:K,
  slope_mean = state_means[, "slope_10_1"],
  infl_mean  = state_means[, "hicp_yoy"],
  dd_mean    = state_means[, "eq_drawdown_12m"],
  ry_mean    = state_means[, "real_yield10"]
)

slope_med_states <- median(state_score$slope_mean)
infl_med_states  <- median(state_score$infl_mean)

state_score$quadrant <- with(state_score,
                             ifelse(slope_mean >= slope_med_states & infl_mean <  infl_med_states, "Normal_expansion",
                                    ifelse(slope_mean >= slope_med_states & infl_mean >= infl_med_states, "Inflationary_exp",
                                           ifelse(slope_mean <  slope_med_states & infl_mean <  infl_med_states, "Late_cycle",
                                                  "Stagflation_risk")))
)

if (anyDuplicated(state_score$quadrant) || any(is.na(state_score$quadrant))) {
  cat("Note: quadrant assignment ambiguous; using fallback rank-based mapping.\n")
  state_score$rank_slope <- rank(-state_score$slope_mean)
  state_score$rank_infl  <- rank(state_score$infl_mean)
  state_score$combined <- state_score$rank_slope + state_score$rank_infl
  ordering <- order(state_score$combined)
  state_score$quadrant <- c("Normal_expansion", "Late_cycle",
                            "Inflationary_exp", "Stagflation_risk")[ordering]
}

cat("\n--- STATE TO REGIME MAPPING ---\n")
print(state_score)

state_to_regime <- setNames(state_score$quadrant, state_score$state)

alloc_map <- c(Normal_expansion = 0.70, Inflationary_exp = 0.55,
               Late_cycle = 0.45, Stagflation_risk = 0.30)

state_to_alloc <- alloc_map[state_to_regime]
names(state_to_alloc) <- 1:K
cat("\n--- STATE TO EQUITY WEIGHT MAPPING ---\n")
print(state_to_alloc)

# ---- 8. Historical Viterbi and posterior probabilities ---------------
cat("\nComputing Viterbi (most likely state sequence) on historical sample...\n")
post_obj <- depmixS4::posterior(fit_main, type = "viterbi")
vit_states <- post_obj$state
cat(sprintf("Viterbi state distribution:\n"))
print(table(vit_states))
cat(sprintf("Distribution shares (%%):\n"))
print(round(prop.table(table(vit_states)) * 100, 1))

post_smooth <- depmixS4::posterior(fit_main, type = "smoothing")
post_probs <- as.matrix(post_smooth)
if (ncol(post_probs) == K + 1) post_probs <- post_probs[, 2:(K + 1), drop = FALSE]
colnames(post_probs) <- paste0("p_state_", 1:K)

row_sums <- rowSums(post_probs)
cat(sprintf("Posterior row sums: min=%.4f max=%.4f (should be 1)\n",
            min(row_sums), max(row_sums)))

hist_alloc_hard <- state_to_alloc[as.character(vit_states)]
hist_alloc_soft <- as.numeric(post_probs %*% as.numeric(state_to_alloc))

cat(sprintf("\nHistorical allocation summary:\n"))
cat(sprintf("  Hard: mean=%.3f, sd=%.3f, range=[%.2f, %.2f]\n",
            mean(hist_alloc_hard), sd(hist_alloc_hard),
            min(hist_alloc_hard), max(hist_alloc_hard)))
cat(sprintf("  Soft: mean=%.3f, sd=%.3f, range=[%.2f, %.2f]\n\n",
            mean(hist_alloc_soft), sd(hist_alloc_soft),
            min(hist_alloc_soft), max(hist_alloc_soft)))

# ---- 9. HMM vs V3 comparison ----------------------------------------
cat("--- HMM vs V3 SUPERVISED REGIME COMPARISON ---\n")
v3_data <- readRDS("Data/features_with_regimes.rds")
v3_regimes <- v3_data$regime
hmm_regimes <- factor(state_to_regime[as.character(vit_states)],
                      levels = c("Normal_expansion", "Inflationary_exp",
                                 "Late_cycle", "Stagflation_risk"))

cat("V3 distribution (counts):\n")
print(table(v3_regimes))
cat("\nHMM distribution (counts):\n")
print(table(hmm_regimes))

cat("\nCross-tabulation (rows = V3, cols = HMM):\n")
print(table(V3 = v3_regimes, HMM = hmm_regimes))

agreement <- mean(as.character(v3_regimes) == as.character(hmm_regimes))
cat(sprintf("\nAgreement V3 vs HMM: %.1f%%\n", 100 * agreement))

# ---- 10. DNB Viterbi+FB with custom decoder (uses RAW scale) ----------
cat("\nBuilding DNB feature trajectories for HMM application...\n")
dnb_features <- array(NA_real_, dim = c(N_SCEN, N_YEARS, 4))
dimnames(dnb_features)[[3]] <- clust_features

for (t in 1:N_YEARS) {
  y1_t  <- y1[,  t + 1]
  y10_t <- y10[, t + 1]
  hicp_t <- hicp_yoy_dnb[, t]
  eq_ret_t <- eq_ret[, t]
  
  slope_t        <- y10_t - y1_t
  real_yield10_t <- y10_t - hicp_t
  eq_dd_t        <- pmin(0, eq_ret_t)
  
  dnb_features[, t, "slope_10_1"]      <- slope_t
  dnb_features[, t, "hicp_yoy"]        <- hicp_t
  dnb_features[, t, "eq_drawdown_12m"] <- eq_dd_t
  dnb_features[, t, "real_yield10"]    <- real_yield10_t
}

cat("DNB feature summary at year 5 (sanity, RAW scale):\n")
for (f in clust_features) {
  cat(sprintf("  %-20s  mean=%.3f  range=[%.3f, %.3f]\n",
              f, mean(dnb_features[, 5, f]),
              min(dnb_features[, 5, f]),
              max(dnb_features[, 5, f])))
}
cat("\n")

cat("Running Viterbi decoding on 20,000 DNB scenarios...\n")
cat("(This may take a few minutes)\n")
start_time <- Sys.time()

init_probs <- fit_main@prior@parameters$coefficients
init_probs <- init_probs / sum(init_probs)

# Custom Viterbi/FB use the RAW-scale state_means and state_sds, which match
# the RAW-scale DNB features. Decoder is identical to V1 (math is correct).
log_emission_prob <- function(x_vec, mean_vec, sd_vec) {
  sum(dnorm(x_vec, mean = mean_vec, sd = sd_vec, log = TRUE))
}

viterbi_decode <- function(obs_seq, trans_mat, state_means, state_sds, init_probs) {
  Tn <- nrow(obs_seq)
  K  <- nrow(state_means)
  log_delta <- matrix(-Inf, Tn, K)
  psi       <- matrix(0,    Tn, K)
  
  log_trans <- log(trans_mat + 1e-12)
  log_init  <- log(init_probs + 1e-12)
  
  for (k in 1:K) {
    log_delta[1, k] <- log_init[k] +
      log_emission_prob(obs_seq[1, ], state_means[k, ], state_sds[k, ])
  }
  for (t in 2:Tn) {
    for (k in 1:K) {
      vec <- log_delta[t - 1, ] + log_trans[, k]
      log_delta[t, k] <- max(vec) +
        log_emission_prob(obs_seq[t, ], state_means[k, ], state_sds[k, ])
      psi[t, k] <- which.max(vec)
    }
  }
  states <- integer(Tn)
  states[Tn] <- which.max(log_delta[Tn, ])
  for (t in (Tn - 1):1) {
    states[t] <- psi[t + 1, states[t + 1]]
  }
  states
}

forward_backward_post <- function(obs_seq, trans_mat, state_means, state_sds, init_probs) {
  Tn <- nrow(obs_seq)
  K  <- nrow(state_means)
  log_alpha <- matrix(-Inf, Tn, K)
  log_beta  <- matrix(-Inf, Tn, K)
  log_trans <- log(trans_mat + 1e-12)
  log_init  <- log(init_probs + 1e-12)
  
  log_em <- matrix(0, Tn, K)
  for (t in 1:Tn) {
    for (k in 1:K) {
      log_em[t, k] <- log_emission_prob(obs_seq[t, ], state_means[k, ], state_sds[k, ])
    }
  }
  log_alpha[1, ] <- log_init + log_em[1, ]
  for (t in 2:Tn) {
    for (k in 1:K) {
      log_alpha[t, k] <- log_em[t, k] +
        matrixStats::logSumExp(log_alpha[t - 1, ] + log_trans[, k])
    }
  }
  log_beta[Tn, ] <- 0
  for (t in (Tn - 1):1) {
    for (k in 1:K) {
      log_beta[t, k] <- matrixStats::logSumExp(
        log_trans[k, ] + log_em[t + 1, ] + log_beta[t + 1, ]
      )
    }
  }
  log_gamma <- log_alpha + log_beta
  for (t in 1:Tn) {
    log_gamma[t, ] <- log_gamma[t, ] - matrixStats::logSumExp(log_gamma[t, ])
  }
  exp(log_gamma)
}

hard_states_dnb <- matrix(NA_integer_, N_SCEN, N_YEARS)
soft_alloc_dnb  <- matrix(NA_real_,    N_SCEN, N_YEARS)

for (i in 1:N_SCEN) {
  obs_seq <- dnb_features[i, , ]
  hard_states_dnb[i, ] <- viterbi_decode(obs_seq, trans_mat,
                                         state_means, state_sds, init_probs)
  post <- forward_backward_post(obs_seq, trans_mat,
                                state_means, state_sds, init_probs)
  soft_alloc_dnb[i, ] <- as.numeric(post %*% as.numeric(state_to_alloc))
  
  if (i %% 2000 == 0) {
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    cat(sprintf("  Scenario %d / %d  (%.1fs elapsed)\n", i, N_SCEN, elapsed))
  }
}

elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
cat(sprintf("\nViterbi+FB decoding completed in %.1fs\n\n", elapsed))

# ---- 11. Build hard/soft DNB equity weight matrices ---------------
hard_alloc_dnb <- matrix(state_to_alloc[as.character(hard_states_dnb)],
                         N_SCEN, N_YEARS)

cat("--- DNB ALLOCATION STATISTICS ---\n")
cat(sprintf("HMM_hard:\n"))
cat(sprintf("  mean=%.3f, sd=%.3f, range=[%.2f, %.2f]\n",
            mean(hard_alloc_dnb), sd(hard_alloc_dnb),
            min(hard_alloc_dnb), max(hard_alloc_dnb)))
cat(sprintf("  Unique values: %s\n",
            paste(sort(unique(as.vector(hard_alloc_dnb))), collapse = ", ")))
cat(sprintf("HMM_soft:\n"))
cat(sprintf("  mean=%.3f, sd=%.3f, range=[%.2f, %.2f]\n\n",
            mean(soft_alloc_dnb), sd(soft_alloc_dnb),
            min(soft_alloc_dnb), max(soft_alloc_dnb)))

cat("HMM_hard DNB state distribution:\n")
print(round(prop.table(table(hard_states_dnb)) * 100, 1))

# ---- 12. Wealth, RR, MDD, turnover ----------------------------------
compute_wealth_path <- function(w_eq) {
  wp <- matrix(0, nrow = N_SCEN, ncol = N_YEARS + 1)
  for (t in 1:N_YEARS) {
    contrib  <- CONTRIB_RATE * salary[, t]
    port_ret <- w_eq[, t] * eq_ret[, t] + (1 - w_eq[, t]) * bd_ret[, t]
    wp[, t + 1] <- (wp[, t] + contrib) * (1 + port_ret)
  }
  wp
}

compute_mdd <- function(wp_matrix) {
  N <- nrow(wp_matrix)
  out <- numeric(N)
  for (i in 1:N) {
    path <- wp_matrix[i, ]
    rm <- cummax(path)
    dd <- ifelse(rm > 0, (path - rm) / rm, 0)
    out[i] <- min(dd)
  }
  out
}

compute_turnover <- function(w_eq) {
  Tn <- ncol(w_eq)
  if (Tn < 2) return(rep(0, nrow(w_eq)))
  rowSums(abs(w_eq[, 2:Tn] - w_eq[, 1:(Tn - 1)]))
}

cat("\nComputing wealth, MDD, turnover, RR for HMM_hard and HMM_soft...\n")
wp_hard <- compute_wealth_path(hard_alloc_dnb)
wp_soft <- compute_wealth_path(soft_alloc_dnb)
mdd_hard <- compute_mdd(wp_hard)
mdd_soft <- compute_mdd(wp_soft)
to_hard <- compute_turnover(hard_alloc_dnb)
to_soft <- compute_turnover(soft_alloc_dnb)
RR_hard <- wp_hard[, N_YEARS + 1] / af / salary_at_retirement
RR_soft <- wp_soft[, N_YEARS + 1] / af / salary_at_retirement

cat(sprintf("HMM_hard:  RR_median=%.3f, shortfall_70=%.3f, CVaR_5=%.3f\n",
            median(RR_hard), mean(RR_hard < 0.70),
            mean(RR_hard[RR_hard <= quantile(RR_hard, 0.05)])))
cat(sprintf("HMM_soft:  RR_median=%.3f, shortfall_70=%.3f, CVaR_5=%.3f\n\n",
            median(RR_soft), mean(RR_soft < 0.70),
            mean(RR_soft[RR_soft <= quantile(RR_soft, 0.05)])))

# ---- 13. Build metrics rows and append --------------------------------
build_metrics_row <- function(strategy_name, RR_vec, mdd_vec, to_vec) {
  data.frame(
    strategy         = strategy_name,
    RR_mean          = round(mean(RR_vec), 3),
    RR_median        = round(median(RR_vec), 3),
    RR_sd            = round(sd(RR_vec), 3),
    RR_p10           = round(quantile(RR_vec, 0.10), 3),
    RR_p25           = round(quantile(RR_vec, 0.25), 3),
    RR_p75           = round(quantile(RR_vec, 0.75), 3),
    RR_p90           = round(quantile(RR_vec, 0.90), 3),
    RR_iqr           = round(quantile(RR_vec, 0.75) - quantile(RR_vec, 0.25), 3),
    P_shortfall_50pct = round(mean(RR_vec < 0.50), 3),
    P_shortfall_60pct = round(mean(RR_vec < 0.60), 3),
    P_shortfall_70pct = round(mean(RR_vec < 0.70), 3),
    P_shortfall_80pct = round(mean(RR_vec < 0.80), 3),
    CVaR_5pct        = round(mean(RR_vec[RR_vec <= quantile(RR_vec, 0.05)]), 3),
    CVaR_10pct       = round(mean(RR_vec[RR_vec <= quantile(RR_vec, 0.10)]), 3),
    RR_sharpe_like   = round(mean(RR_vec) / sd(RR_vec), 3),
    MDD_mean         = round(mean(mdd_vec), 3),
    MDD_median       = round(median(mdd_vec), 3),
    MDD_p10          = round(quantile(mdd_vec, 0.10), 3),
    turnover_mean    = round(mean(to_vec), 3),
    turnover_median  = round(median(to_vec), 3),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

row_hard <- build_metrics_row("HMM_hard", RR_hard, mdd_hard, to_hard)
row_soft <- build_metrics_row("HMM_soft", RR_soft, mdd_soft, to_soft)

metrics_new <- bind_rows(metrics_old, row_hard, row_soft)

wp_new <- wp_old
wp_new[["HMM_hard"]] <- wp_hard
wp_new[["HMM_soft"]] <- wp_soft

mdd_new <- cbind(mdd_old, HMM_hard = mdd_hard, HMM_soft = mdd_soft)
to_new  <- cbind(to_old,  HMM_hard = to_hard,  HMM_soft = to_soft)
RR_new  <- cbind(RR_old,  HMM_hard = RR_hard,  HMM_soft = RR_soft)

# ---- 14. Comparison display ----------------------------------------
cat("--- METRICS COMPARISON ---\n")
key_strategies <- c("60/40", "Glide", "MV_plain", "RF",
                    "HybridNaive", "Direct_PPP", "HMM_hard", "HMM_soft")
print(metrics_new %>%
        filter(strategy %in% key_strategies) %>%
        select(strategy, RR_median, RR_sd, P_shortfall_70pct,
               CVaR_5pct, MDD_mean, RR_sharpe_like, turnover_mean),
      row.names = FALSE)

# ---- 15. Save outputs ----------------------------------------------
saveRDS(list(
  fit = fit_main,
  state_means_raw = state_means,
  state_sds_raw = state_sds,
  state_means_std = state_means_std,
  state_sds_std = state_sds_std,
  feat_mean = feat_mean,
  feat_sd = feat_sd,
  trans_mat = trans_mat,
  init_probs = init_probs,
  state_to_regime = state_to_regime,
  state_to_alloc = state_to_alloc,
  bic_diagnostic = bic_results,
  K_selected = K,
  K_best_bic = best_k_bic,
  hist_viterbi_states = vit_states,
  hist_posterior_probs = post_probs,
  v3_hmm_agreement = agreement
), "Data/hmm_model.rds")

saveRDS(list(
  hard = hist_alloc_hard,
  soft = hist_alloc_soft,
  v3   = as.character(v3_regimes)
), "Data/hmm_state_assignments.rds")

saveRDS(list(
  hard = hard_alloc_dnb,
  soft = soft_alloc_dnb,
  states = hard_states_dnb
), "Data/hmm_dnb_allocations.rds")

saveRDS(metrics_new, "Data/metrics_full.rds")
saveRDS(wp_new,      "Data/wealth_paths.rds")
saveRDS(mdd_new,     "Data/max_drawdowns.rds")
saveRDS(to_new,      "Data/turnovers.rds")
saveRDS(RR_new,      "Data/replacement_ratios_with_hybrid_naive.rds")

cat("\nSaved (HMM unsupervised regimes, V2 standardized):\n")
cat("  Data/hmm_model.rds                  (fit + diagnostics)\n")
cat("  Data/hmm_state_assignments.rds      (historical labels)\n")
cat("  Data/hmm_dnb_allocations.rds        (20000 x 42 x 2)\n")
cat("  Data/metrics_full.rds               (12 strategies)\n")
cat("  Data/wealth_paths.rds               (12 strategies)\n")
cat("  Data/max_drawdowns.rds              (12 strategies)\n")
cat("  Data/turnovers.rds                  (12 strategies)\n")
cat("  Data/replacement_ratios_with_hybrid_naive.rds (12 strategies)\n")
