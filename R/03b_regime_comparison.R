# 03b_regime_comparison.R
#
# Input:        Data/features_historical.rds
# Output:       Console output only (no save until we pick a variant)
# Dependencies: dplyr, tidyr, lubridate, zoo
#
# Compares 3 regime taxonomy variants:
#   V1 - Contemporaneous slope and inflation (current)
#   V2 - Lagged slope (12 months) and contemporaneous inflation
#   V3 - V1 with equity drawdown override

library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)

df <- readRDS("Data/features_historical.rds")

INFL_THRESHOLD  <- 2.0
SLOPE_THRESHOLD <- 1.0
DRAWDOWN_THRESHOLD <- -0.15  # equity drawdown threshold for V3

# ---- Compute rolling 12-month equity drawdown for V3 ----------------------
df <- df %>%
  arrange(date) %>%
  mutate(
    eq_max_12m = rollapply(eq_price, width = 12, FUN = max, align = "right", fill = NA),
    eq_drawdown = (eq_price / eq_max_12m) - 1
  )

# ---- Helper to assign regime given slope and inflation --------------------
assign_regime <- function(slope, infl) {
  case_when(
    slope >  SLOPE_THRESHOLD & infl <  INFL_THRESHOLD ~ "Normal_expansion",
    slope >  SLOPE_THRESHOLD & infl >= INFL_THRESHOLD ~ "Inflationary_exp",
    slope <= SLOPE_THRESHOLD & infl <  INFL_THRESHOLD ~ "Late_cycle",
    slope <= SLOPE_THRESHOLD & infl >= INFL_THRESHOLD ~ "Stagflation_risk"
  )
}

# ---- Build the 3 variants -------------------------------------------------
df <- df %>%
  mutate(
    regime_v1 = assign_regime(slope_10_1, hicp_yoy),                 # contemporaneous
    regime_v2 = assign_regime(slope_l12,  hicp_yoy),                 # lagged slope
    regime_v3 = case_when(
      !is.na(eq_drawdown) & eq_drawdown < DRAWDOWN_THRESHOLD ~ "Late_cycle",
      TRUE ~ regime_v1
    )
  ) %>%
  mutate(across(c(regime_v1, regime_v2, regime_v3),
                ~ factor(., levels = c("Normal_expansion", "Inflationary_exp",
                                       "Late_cycle", "Stagflation_risk"))))

# ---- V2 needs slope_l12: keep only rows where it exists --------------------
# (slope_l12 was already in the feature set so this is just a check)
df_v <- df %>% filter(!is.na(regime_v2))
cat("Observations available for comparison:", nrow(df_v), "\n\n")

# ---- 1. Class distribution ------------------------------------------------
cat("=== CLASS DISTRIBUTION ===\n\n")
for (v in c("regime_v1", "regime_v2", "regime_v3")) {
  cat(toupper(v), ":\n")
  print(table(df_v[[v]]))
  cat("\nShares (%):\n")
  print(round(prop.table(table(df_v[[v]])) * 100, 1))
  cat("\n---\n")
}

# ---- 2. Key dates comparison ----------------------------------------------
cat("\n=== KEY DATES ===\n")
key_dates <- as.Date(c(
  "2007-12-31",  # pre-crisis, curve flattening
  "2008-12-31",  # Lehman crisis
  "2009-06-30",  # crisis recovery
  "2011-12-31",  # debt crisis
  "2013-06-30",  # recovery
  "2020-06-30",  # COVID
  "2022-06-30",  # inflation surge
  "2022-12-31",  # peak inflation
  "2023-12-31"   # disinflation
))

cmp <- df_v %>%
  filter(date %in% key_dates) %>%
  select(date, hicp_yoy, slope_10_1, slope_l12, eq_drawdown,
         regime_v1, regime_v2, regime_v3)
print(cmp, n = 30)

# ---- 3. Transition rate ---------------------------------------------------
cat("\n=== TRANSITION RATES ===\n")
for (v in c("regime_v1", "regime_v2", "regime_v3")) {
  reg <- df_v[[v]]
  reg_prev <- lag(reg)
  n_trans <- sum(reg != reg_prev, na.rm = TRUE)
  n_obs   <- sum(!is.na(reg_prev))
  cat(sprintf("%-12s: %d / %d (%.1f%%)\n",
              v, n_trans, n_obs, 100 * n_trans / n_obs))
}

# ---- 4. Allocation during known crisis windows ----------------------------
allocation_map <- c(
  Normal_expansion  = 0.70,
  Inflationary_exp  = 0.55,
  Late_cycle        = 0.45,
  Stagflation_risk  = 0.30
)

crisis_windows <- list(
  "GFC 2008-2009"    = c("2008-01-01", "2009-12-31"),
  "Debt crisis 2011" = c("2011-01-01", "2012-12-31"),
  "COVID 2020"       = c("2020-01-01", "2020-12-31"),
  "Inflation 2022"   = c("2022-01-01", "2023-06-30")
)

cat("\n=== MEAN EQUITY ALLOCATION DURING CRISIS WINDOWS ===\n")
for (cn in names(crisis_windows)) {
  win <- as.Date(crisis_windows[[cn]])
  sub <- df_v %>% filter(date >= win[1] & date <= win[2])
  if (nrow(sub) == 0) next
  alloc <- sapply(c("regime_v1", "regime_v2", "regime_v3"), function(v) {
    mean(allocation_map[as.character(sub[[v]])], na.rm = TRUE)
  })
  cat(sprintf("%-20s (n=%2d):  V1=%.2f  V2=%.2f  V3=%.2f\n",
              cn, nrow(sub), alloc[1], alloc[2], alloc[3]))
}
