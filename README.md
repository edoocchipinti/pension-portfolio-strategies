# Robust Long-Term Portfolio Strategies for Pension Funds

R code for my MSc thesis in Financial Economics (Maastricht University, School of
Business and Economics, 2026).

The project compares **model-based**, **data-driven** and **model-free** asset
allocation strategies for a long-horizon pension fund investor, and evaluates
how robust each one is to estimation error and macroeconomic uncertainty.

## Setup

A single participant is followed over a **42-year accumulation phase** across
**20,000 DNB scenarios**, contributing a fixed share of salary each year into a
two-asset portfolio (global equity / 10-year bonds). The outcome variable is the
**replacement ratio (RR)** at retirement: accumulated wealth converted into a
lifetime annuity using a scenario-specific annuity factor built from AG2024
cohort survival probabilities and the scenario's own yield curve.

Strategies are therefore judged on the pension outcome they deliver, not on
portfolio returns in isolation.

## Data

- **DNB uniform scenario set (CP2022)** — 20,000 scenarios of equity returns,
  inflation and the yield curve (1y, 5y, 10y, 20y).
- **Historical monthly series (EU)** — MSCI World total return, HICP, Euribor 3m,
  Bund yields at 1y/5y/10y/20y; plus VSTOXX, M3 and unemployment, retained for
  robustness checks but excluded from the feature set because they are not
  reproducible inside the DNB scenarios.
- **AG2024 cohort survival table** for the annuity factor.

Raw data is not redistributed here. See [`data/README.md`](data/README.md) for
sources and the expected file layout.

## Strategies compared

| Family | Strategies |
|---|---|
| Fixed rules | `60/40`, `1/N`, `Glide` (80% → 30% equity) |
| Model-based | `MV_plain` (mean-variance), `MV_LW` (Ledoit-Wolf shrinkage) |
| Data-driven (ML) | `Lasso`, `Ridge`, `RF` (random forest) on 23 macro-financial features |
| Regime-based | Supervised rule labels; unsupervised Gaussian HMM (`HMM_hard` via Viterbi, `HMM_soft` via forward-backward) |
| End-to-end | `Direct_PPP` — parametric portfolio policy, Brandt, Santa-Clara & Valkanov (2009) |
| Combinations | `Hybrid`, `HybridNaive` |

Two design choices matter for the comparison to be honest:

- **No look-ahead.** ML allocations use out-of-fold predictions and lagged
  features only; the look-ahead present in the first version of the hybrid and
  direct-policy scripts is documented and corrected in `11_` and `12_`.
- **Regime count chosen ex ante.** The HMM moves from K=4 to K=3 on BIC stability
  and EM convergence grounds, declared before comparing performance, not after.

## Evaluation

For each strategy: median RR, dispersion (sd, IQR, p10–p90), shortfall
probability against 50/60/70/80% targets, CVaR at 5% and 10%, maximum drawdown
of the wealth path, and turnover.

## Robustness

Each `10*` script perturbs one assumption and re-runs the full pipeline:
allocation map, regime thresholds, drawdown trigger, bond duration, mean-variance
estimation window, risk aversion, contribution rate, equity cap, wage drift.

## Repository structure

| Script | Purpose |
|---|---|
| `01_compute_annuity_factor.R` | Scenario-specific annuity factor from AG2024 survival + scenario yield curve |
| `02_feature_engineering.R` | Builds the 23-feature panel from historical monthly data |
| `03_regime_labels.R` / `03b_` | Supervised regime labels and comparison of labelling rules |
| `04_train_ml.R` | Trains Lasso, Ridge and random forest with cross-validation |
| `05_calibrate_mv.R` / `05b_` | Mean-variance calibration; DNB yield extraction |
| `06_apply_to_dnb.R` | Maps historical models onto the 20,000 DNB scenarios |
| `07_compute_rr.R` | Wealth accumulation and replacement ratios |
| `08_metrics.R` | Master metrics table across strategies |
| `09_plots.R` / `09b_` | Figures; appends hybrid strategies to the metrics table |
| `10a`–`10q` | Robustness checks (one assumption each) |
| `11_hybrid_strategy.R` | Hybrid strategy on out-of-fold ML predictions |
| `12_direct_ml.R` / `12b_` | Direct parametric portfolio policy |
| `13_unsupervised_regimes_k4.R` | HMM with K=4 (superseded, kept for the robustness discussion) |
| `13b_unsupervised_regimes_k3.R` / `13b_hmm_recompute.R` | Final HMM specification with K=3 |
| `15_export_results.R` | Exports the result tables used in the thesis |
| `run_all.R` | Runs the pipeline end to end |

## Requirements

R 4.x and the following packages:

```r
install.packages(c("readxl", "dplyr", "tidyr", "lubridate", "zoo",
                   "glmnet", "ranger", "depmixS4", "matrixStats",
                   "ggplot2", "scales", "patchwork", "RColorBrewer"))
```

## Author

Edoardo Occhipinti — MSc Financial Economics, Maastricht University
[linkedin.com/in/edo-occhipinti](https://linkedin.com/in/edo-occhipinti)
