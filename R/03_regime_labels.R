# 03_regime_labels.R  (V3 — drawdown override)
#
# Input:        Data/features_historical.rds
# Output:       Data/features_with_regimes.rds
#               Data/allocation_map.rds
# Dependencies: dplyr, tidyr, lubridate
#
# Assigns 4 macro-financial regimes using a two-stage logic:
#
# Stage 1 (contemporaneous slope x inflation):
#   Normal expansion       : slope >  1.0%, HICP YoY <  2.0%
#   Inflationary expansion : slope >  1.0%, HICP YoY >= 2.0%
#   Late cycle             : slope <= 1.0%, HICP YoY <  2.0%
#   Stagflation risk       : slope <= 1.0%, HICP YoY >= 2.0%
#
# Stage 2 (drawdown override):
#   If eq_drawdown_12m < -15%, regime is forced to Late_cycle.
#   This corrects the failure of Stage 1 in post-crisis recovery periods
#   (2008-12, 2009-06) where the yield curve steepens due to emergency
#   rate cuts while the real economy is still contracting.
#
# Allocation mapping (pre-fixed, not estimated):
#   Normal expansion       -> 70% equity
#   Inflationary expansion -> 55% equity
#   Late cycle             -> 45% equity
#   Stagflation risk       -> 30% equity
#
# Rationale and citations:
#   Inflation threshold:   ECB price stability target (2%)
#   Slope threshold:       NY Fed convention (Adrian, Crump & Moench, 2010)
#   Drawdown override:     bear-market correction signal (Faber 2007;
#                          Estrada 2008)
#   Drawdown threshold:    15% — conventional "bear correction" magnitude

library(dplyr)
library(tidyr)
library(lubridate)

# ---- 1. Load features -----------------------------------------------------
df <- readRDS("Data/features_historical.rds")

# ---- 2. Thresholds --------------------------------------------------------
INFL_THRESHOLD     <- 2.0
SLOPE_THRESHOLD    <- 1.0
DRAWDOWN_THRESHOLD <- -0.15

# ---- 3. Stage 1: base regime from slope x inflation -----------------------
df <- df %>%
  mutate(
    regime_base = case_when(
      slope_10_1 >  SLOPE_THRESHOLD & hicp_yoy <  INFL_THRESHOLD ~ "Normal_expansion",
      slope_10_1 >  SLOPE_THRESHOLD & hicp_yoy >= INFL_THRESHOLD ~ "Inflationary_exp",
      slope_10_1 <= SLOPE_THRESHOLD & hicp_yoy <  INFL_THRESHOLD ~ "Late_cycle",
      slope_10_1 <= SLOPE_THRESHOLD & hicp_yoy >= INFL_THRESHOLD ~ "Stagflation_risk"
    )
  )

# ---- 4. Stage 2: drawdown override ----------------------------------------
df <- df %>%
  mutate(
    drawdown_override_active = eq_drawdown_12m < DRAWDOWN_THRESHOLD,
    regime = ifelse(drawdown_override_active, "Late_cycle", regime_base),
    regime = factor(regime,
                    levels = c("Normal_expansion", "Inflationary_exp",
                               "Late_cycle", "Stagflation_risk")),
    regime_base = factor(regime_base,
                         levels = c("Normal_expansion", "Inflationary_exp",
                                    "Late_cycle", "Stagflation_risk"))
  )

# ---- 5. Allocation map ----------------------------------------------------
allocation_map <- c(
  Normal_expansion  = 0.70,
  Inflationary_exp  = 0.55,
  Late_cycle        = 0.45,
  Stagflation_risk  = 0.30
)

df <- df %>%
  mutate(target_equity = allocation_map[as.character(regime)])

# ---- 6. Override impact ---------------------------------------------------
cat("--- OVERRIDE IMPACT ---\n")
n_override <- sum(df$drawdown_override_active)
cat(sprintf("Observations with active override: %d / %d (%.1f%%)\n",
            n_override, nrow(df), 100 * n_override / nrow(df)))

cat("\nReclassifications by base regime (base -> Late_cycle):\n")
print(df %>%
        filter(drawdown_override_active) %>%
        count(regime_base, name = "n"))

# ---- 7. Class distribution: base vs final --------------------------------
cat("\n--- BASE (V1) VS FINAL (V3) DISTRIBUTION ---\n")
cat("\nBase regime (V1):\n")
print(table(df$regime_base))
cat("\nFinal regime (V3):\n")
print(table(df$regime))
cat("\nFinal shares (%):\n")
print(round(prop.table(table(df$regime)) * 100, 1))

class_counts <- as.numeric(table(df$regime))
if (any(class_counts < 30)) {
  cat("\nWARNING: at least one regime has < 30 observations.\n")
} else {
  cat("\nAll regimes have >= 30 observations. Class balance OK.\n")
}

# ---- 8. Regime distribution by year --------------------------------------
cat("\n--- REGIME BY YEAR (V3 final) ---\n")
df_year <- df %>%
  mutate(year = year(date)) %>%
  count(year, regime) %>%
  pivot_wider(names_from = regime, values_from = n, values_fill = 0) %>%
  arrange(year)
print(df_year, n = 30)

# ---- 9. Sanity check on key historical dates -----------------------------
cat("\n--- SANITY CHECK: KEY DATES ---\n")
sanity_dates <- as.Date(c(
  "2007-12-31", "2008-12-31", "2009-06-30",
  "2011-12-31", "2013-06-30", "2020-06-30",
  "2022-06-30", "2022-12-31", "2023-12-31"
))

sanity <- df %>%
  filter(date %in% sanity_dates) %>%
  select(date, hicp_yoy, slope_10_1, eq_drawdown_12m,
         regime_base, regime, target_equity)
print(sanity, n = 20)

# ---- 10. Transition rate -------------------------------------------------
cat("\n--- TRANSITIONS ---\n")
n_trans <- sum(df$regime != lag(df$regime), na.rm = TRUE)
n_obs   <- sum(!is.na(lag(df$regime)))
cat(sprintf("Month-to-month transitions: %d / %d (%.1f%%)\n",
            n_trans, n_obs, 100 * n_trans / n_obs))

# ---- 11. Save -------------------------------------------------------------
saveRDS(df, "Data/features_with_regimes.rds")
saveRDS(allocation_map, "Data/allocation_map.rds")

cat("\nSaved Data/features_with_regimes.rds\n")
cat("Saved Data/allocation_map.rds\n")
