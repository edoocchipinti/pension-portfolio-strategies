# 02_feature_engineering.R  (V3 — adds eq_drawdown_12m, removes IP)
#
# Input:        Data/historical_monthly.xlsx (now extended to Dec 2025)
# Output:       Data/features_historical.rds
#               Data/feature_cols.rds
# Dependencies: readxl, dplyr, lubridate, zoo
#
# Builds 23 features for ML and MV strategies, restricted to variables
# computable both historically and from DNB CP2022 scenarios:
#   - equity returns + drawdown
#   - EU inflation
#   - yield curve (1y, 5y, 10y, 20y)
#
# Industrial Production (IP) is DROPPED entirely because it has a ~45-day
# publication lag which makes the last 3 months of 2025 unavailable, and
# because IP is not a DNB-compatible variable anyway.
#
# Non-DNB-compatible variables retained (VSTOXX, M3, Unemployment) for
# potential robustness checks, but NOT in feature_cols.
#
# All HICP values are assumed available at end-of-month t (flash release
# convention). Yields and equity are end-of-month market data.

library(readxl)
library(dplyr)
library(lubridate)
library(zoo)

# ---- 1. Load and rename columns ------------------------------------------
raw <- read_excel("Data/historical_monthly.xlsx", sheet = 1)

df <- raw %>%
  transmute(
    date       = as.Date(Date),
    eq_price   = `MSCI World_Total Return`,
    hicp       = HICP_EU,
    euribor_3m = `Euribor 3m`,
    y1         = `Bund 1y`,
    y5         = `Bund 5y`,
    y10        = `Bund 10y`,
    y20        = `Bund 20y`,
    # Non-DNB-compatible (kept but not used as ML features)
    vstoxx     = VSTOXX_Price,
    m3         = M3_OBS.VALUE,
    unemp      = `Unemployment rate`
  ) %>%
  arrange(date)

# ---- 2. Base transformations ---------------------------------------------
df <- df %>%
  mutate(
    eq_ret   = log(eq_price / lag(eq_price, 1)),         # log return, monthly
    hicp_yoy = 100 * (log(hicp) - log(lag(hicp, 12))),   # YoY inflation, pp
    dy1      = 100 * (y1  - lag(y1,  1)),                # Δ short yield, bps
    dy10     = 100 * (y10 - lag(y10, 1))                 # Δ long yield, bps
  )

# ---- 3. Yield curve features ---------------------------------------------
df <- df %>%
  mutate(
    slope_10_1   = y10 - y1,            # term spread (Estrella-Mishkin 1998)
    slope_20_5   = y20 - y5,            # long-end slope
    curvature    = 2 * y5 - y1 - y10,   # butterfly
    real_yield10 = y10 - hicp_yoy       # ex-post real long yield
  )

# ---- 4. Equity drawdown --------------------------------------------------
# Rolling 12-month max of equity price; drawdown = current/max - 1.
# Captures cumulative equity stress within the last year.
# Used both as ML feature and (in 03_regime_labels.R) as override trigger.
df <- df %>%
  mutate(
    eq_max_12m       = rollapply(eq_price, width = 12, FUN = max,
                                 align = "right", fill = NA),
    eq_drawdown_12m  = (eq_price / eq_max_12m) - 1
  ) %>%
  select(-eq_max_12m)

# ---- 5. Selected temporal lags -------------------------------------------
df <- df %>%
  mutate(
    eq_ret_l1    = lag(eq_ret, 1),
    eq_ret_l3    = lag(eq_ret, 3),
    eq_ret_l12   = lag(eq_ret, 12),
    hicp_yoy_l3  = lag(hicp_yoy, 3),
    hicp_yoy_l12 = lag(hicp_yoy, 12),
    slope_l3     = lag(slope_10_1, 3),
    slope_l12    = lag(slope_10_1, 12),
    y10_l12      = lag(y10, 12)
  )

# ---- 6. Economically motivated interactions ------------------------------
df <- df %>%
  mutate(
    slope_x_infl = slope_10_1 * hicp_yoy,
    realy_x_eq   = real_yield10 * eq_ret_l1
  )

# ---- 7. Define final feature set and drop NAs ----------------------------
feature_cols <- c(
  # Base (8)
  "eq_ret", "hicp_yoy", "y1", "y5", "y10", "y20", "dy1", "dy10",
  # Yield curve derived (4)
  "slope_10_1", "slope_20_5", "curvature", "real_yield10",
  # Equity drawdown (1)
  "eq_drawdown_12m",
  # Lags (8)
  "eq_ret_l1", "eq_ret_l3", "eq_ret_l12",
  "hicp_yoy_l3", "hicp_yoy_l12",
  "slope_l3", "slope_l12", "y10_l12",
  # Interactions (2)
  "slope_x_infl", "realy_x_eq"
)

cat("Observations before NA drop:", nrow(df), "\n")
df_clean <- df %>% filter(if_all(all_of(feature_cols), ~ !is.na(.)))
cat("Observations after NA drop: ", nrow(df_clean), "\n")
cat("First usable date:          ", as.character(min(df_clean$date)), "\n")
cat("Last usable date:           ", as.character(max(df_clean$date)), "\n")
cat("Number of ML features:      ", length(feature_cols), "\n\n")

# ---- 8. Save -------------------------------------------------------------
saveRDS(df_clean, "Data/features_historical.rds")
saveRDS(feature_cols, "Data/feature_cols.rds")
cat("Saved Data/features_historical.rds\n")
cat("Saved Data/feature_cols.rds\n\n")

# ---- 9. Feature summary --------------------------------------------------
cat("--- FEATURE SUMMARY ---\n")
print(summary(df_clean[, feature_cols]))

cat("\n--- DRAWDOWN STATISTICS ---\n")
cat(sprintf("Min drawdown:        %.2f%%\n", 100 * min(df_clean$eq_drawdown_12m)))
cat(sprintf("Max drawdown:        %.2f%%\n", 100 * max(df_clean$eq_drawdown_12m)))
cat(sprintf("Obs with dd < -15%%:  %d (%.1f%%)\n",
            sum(df_clean$eq_drawdown_12m < -0.15),
            100 * mean(df_clean$eq_drawdown_12m < -0.15)))
