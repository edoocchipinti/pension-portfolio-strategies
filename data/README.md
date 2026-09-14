# Data

Raw data is not included in this repository. The scripts expect the following
files in a `Data/` folder at the project root.

## DNB scenario set (CP2022)

Published by De Nederlandsche Bank as the uniform scenario set for pension
funds. 20,000 scenarios.

| File | Contents |
|---|---|
| `Data/equity_returns.rds` | 20,000 x 42 — annual equity returns |
| `Data/bond_returns_10y.rds` | 20,000 x 42 — annual 10-year bond returns |
| `Data/yield_1y.rds`, `yield_5y.rds`, `yield_10y.rds`, `yield_20y.rds` | 20,000 x 101 — yield curve, decimal |
| `Data/salary_path.rds` | 20,000 x 42 — salary path |

## Historical monthly series

`Data/historical_monthly.xlsx`, one sheet, monthly frequency, with columns:

`Date`, `MSCI World_Total Return`, `HICP_EU`, `Euribor 3m`,
`Bund 1y`, `Bund 5y`, `Bund 10y`, `Bund 20y`,
`VSTOXX_Price`, `M3_OBS.VALUE`, `Unemployment rate`

Sources: MSCI, Eurostat (HICP, unemployment), ECB Statistical Data Warehouse
(Euribor, M3, Bund yields), STOXX (VSTOXX).

## Mortality

`Data/AG2024_cohort_survival.rds` — cohort survival probabilities from the
Koninklijk Actuarieel Genootschap AG2024 projection table.

## Generated files

Everything else under `Data/` is produced by the pipeline (`features_historical.rds`,
`mv_weights.rds`, `dnb_equity_allocations.rds`, `metrics_full.rds`, and so on)
and does not need to be supplied.
