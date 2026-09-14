# run_all.R
#
# Runs the full pipeline in numbered order.
# Expects the raw data files described in data/README.md under Data/.
#
# Note: scripts are ordered by file number. Verify the order on the first
# end-to-end run — some stages depend on outputs written by later-numbered
# scripts if the pipeline is re-run after a partial execution.

scripts <- c(
  "R/05b_dnb_yields.R",
  "R/01_compute_annuity_factor.R",
  "R/02_feature_engineering.R",
  "R/03_regime_labels.R",
  "R/03b_regime_comparison.R",
  "R/04_train_ml.R",
  "R/05_calibrate_mv.R",
  "R/06_apply_to_dnb.R",
  "R/07_compute_rr.R",
  "R/08_metrics.R",
  "R/11_hybrid_strategy.R",
  "R/09b_add_hybrid_to_metrics.R",
  "R/12_direct_ml.R",
  "R/12b_direct_recompute.R",
  "R/13b_unsupervised_regimes_k3.R",
  "R/13b_hmm_recompute.R",
  "R/09_plots.R",
  "R/15_export_results.R"
)

robustness <- c(
  "R/10a_robustness_allocation_map.R",
  "R/10b_robustness_regime_thresholds.R",
  "R/10c_robustness_drawdown_threshold.R",
  "R/10d_robustness_bond_duration.R",
  "R/10e_robustness_mv_window.R",
  "R/10h_robustness_gamma.R",
  "R/10i_robustness_contribution.R",
  "R/10k_robustness_equity_cap.R",
  "R/10q_robustness_wage_drift.R"
)

for (s in c(scripts, robustness)) {
  cat("\n=====", s, "=====\n")
  source(s)
}

cat("\nPipeline complete.\n")
