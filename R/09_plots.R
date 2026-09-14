# 09_plots.R
#
# Input:  Data/metrics_full.rds                  (9 strategies)
#         Data/wealth_paths.rds                  (9 strategies)
#         Data/max_drawdowns.rds                 (9 strategies)
#         Data/replacement_ratios_with_hybrid_naive.rds (20000 x 9)
#         Data/dnb_equity_allocations.rds        (20000 x 42 x 5 for MV/ML)
#         Data/mv_weights.rds
#
# Output: Output/Plots/01_density_facet.pdf
#         Output/Plots/02_shortfall_bar.pdf
#         Output/Plots/03_allocation_evolution.pdf
#         Output/Plots/04_wealth_fan.pdf
#         Output/Plots/05_tradeoff_scatter.pdf
#         Output/Plots/06_hybrid_composition.pdf

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(patchwork)

# ---- 0. Setup ------------------------------------------------------------
if (!dir.exists("Output"))       dir.create("Output")
if (!dir.exists("Output/Plots")) dir.create("Output/Plots")

metrics <- readRDS("Data/metrics_full.rds")
wp      <- readRDS("Data/wealth_paths.rds")
mdd     <- readRDS("Data/max_drawdowns.rds")
RR      <- readRDS("Data/replacement_ratios_with_hybrid_naive.rds")
allocs  <- readRDS("Data/dnb_equity_allocations.rds")
mv_w    <- readRDS("Data/mv_weights.rds")

strategy_order <- c("60/40", "1/N", "Glide", "MV_plain", "MV_LW",
                    "Lasso", "Ridge", "RF", "HybridNaive")

strategy_family <- c(
  "60/40" = "Rule-based", "1/N" = "Rule-based", "Glide" = "Rule-based",
  "MV_plain" = "Statistical", "MV_LW" = "Statistical",
  "Lasso" = "Machine Learning", "Ridge" = "Machine Learning",
  "RF" = "Machine Learning",
  "HybridNaive" = "Hybrid"
)

family_colors <- c(
  "Rule-based"       = "#1f77b4",
  "Statistical"      = "#2ca02c",
  "Machine Learning" = "#d62728",
  "Hybrid"           = "#9467bd"
)

theme_thesis <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.major.y = element_line(color = "grey92"),
      strip.background   = element_rect(fill = "grey95", color = NA),
      strip.text         = element_text(face = "bold"),
      legend.position    = "bottom",
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(color = "grey40", size = 10),
      plot.caption       = element_text(color = "grey50", size = 8, hjust = 0)
    )
}

# ============================================================
# PLOT 1 — Density of RR by strategy (faceted)
# ============================================================
RR_df <- as.data.frame(RR) %>%
  pivot_longer(everything(), names_to = "strategy", values_to = "RR") %>%
  mutate(strategy = factor(strategy, levels = strategy_order),
         family   = strategy_family[as.character(strategy)])

p1 <- ggplot(RR_df, aes(x = RR, fill = family)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_vline(xintercept = 0.70, linetype = "dashed",
             color = "black", linewidth = 0.4) +
  facet_wrap(~ strategy, ncol = 3, scales = "fixed") +
  scale_fill_manual(values = family_colors) +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.5),
                     labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Distribution of Replacement Ratio across 20,000 DNB scenarios",
    subtitle = "Dashed vertical line: pension target (70%)",
    x        = "Replacement Ratio",
    y        = "Density",
    fill     = "Strategy family",
    caption  = "Source: simulation on DNB CP2022 scenarios; 42-year accumulation phase."
  ) +
  theme_thesis()

ggsave("Output/Plots/01_density_facet.pdf", p1,
       width = 9, height = 8, device = "pdf")
cat("Saved 01_density_facet.pdf\n")

# ============================================================
# PLOT 2 — Shortfall bar chart
# ============================================================
shortfall_df <- metrics %>%
  select(strategy, P_shortfall_70pct) %>%
  mutate(
    strategy = factor(strategy, levels = strategy_order),
    family   = strategy_family[as.character(strategy)]
  ) %>%
  arrange(P_shortfall_70pct) %>%
  mutate(strategy_ord = factor(strategy, levels = strategy))

p2 <- ggplot(shortfall_df,
             aes(x = P_shortfall_70pct, y = strategy_ord, fill = family)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = percent(P_shortfall_70pct, accuracy = 0.1)),
            hjust = -0.1, size = 3.3) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = family_colors) +
  labs(
    title    = "Probability of shortfall below 70% replacement ratio",
    subtitle = "Lower is better. Strategies ranked from most to least robust.",
    x        = "P(RR < 0.70)",
    y        = NULL,
    fill     = "Strategy family",
    caption  = "Source: simulation on DNB CP2022 scenarios."
  ) +
  theme_thesis() +
  theme(panel.grid.major.y = element_blank())

ggsave("Output/Plots/02_shortfall_bar.pdf", p2,
       width = 8, height = 5, device = "pdf")
cat("Saved 02_shortfall_bar.pdf\n")

# ============================================================
# PLOT 3 — Allocation evolution
# ============================================================
N_SCEN  <- dim(allocs)[1]
N_YEARS <- dim(allocs)[2]
year_seq <- 1:N_YEARS

w_6040  <- matrix(0.60, N_SCEN, N_YEARS)
w_1N    <- matrix(0.50, N_SCEN, N_YEARS)
w_glide <- matrix(rep(seq(0.80, 0.30, length.out = N_YEARS),
                      each = N_SCEN), N_SCEN, N_YEARS)
w_mv_p  <- matrix(as.numeric(mv_w$plain["equity"]), N_SCEN, N_YEARS)
w_mv_lw <- matrix(as.numeric(mv_w$lw["equity"]),    N_SCEN, N_YEARS)
w_lasso <- allocs[, , "Lasso"]
w_ridge <- allocs[, , "Ridge"]
w_rf    <- allocs[, , "RF"]
w_hybrid <- (1/3) * w_6040 + (1/3) * w_mv_p + (1/3) * w_rf

w_list <- list(
  `60/40` = w_6040, `1/N` = w_1N, Glide = w_glide,
  MV_plain = w_mv_p, MV_LW = w_mv_lw,
  Lasso = w_lasso, Ridge = w_ridge, RF = w_rf,
  HybridNaive = w_hybrid
)

alloc_evo <- lapply(names(w_list), function(s) {
  data.frame(
    strategy = s, year = year_seq,
    mean_w   = colMeans(w_list[[s]]),
    p25_w    = apply(w_list[[s]], 2, quantile, 0.25),
    p75_w    = apply(w_list[[s]], 2, quantile, 0.75)
  )
}) %>% bind_rows() %>%
  mutate(strategy = factor(strategy, levels = strategy_order),
         family   = strategy_family[as.character(strategy)])

p3 <- ggplot(alloc_evo, aes(x = year, y = mean_w, color = family)) +
  geom_ribbon(aes(ymin = p25_w, ymax = p75_w, fill = family),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ strategy, ncol = 3) +
  scale_color_manual(values = family_colors) +
  scale_fill_manual(values = family_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.25, 0.85), breaks = seq(0.3, 0.8, 0.1)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 42)) +
  labs(
    title    = "Equity allocation over the 42-year accumulation horizon",
    subtitle = "Line: cross-scenario mean. Shaded band: 25th-75th percentile.",
    x        = "Career year (age 25 = year 1)",
    y        = "Equity weight",
    color    = "Strategy family",
    fill     = "Strategy family",
    caption  = "Static strategies show flat lines. ML strategies show scenario-dependent variability."
  ) +
  theme_thesis()

ggsave("Output/Plots/03_allocation_evolution.pdf", p3,
       width = 9, height = 8, device = "pdf")
cat("Saved 03_allocation_evolution.pdf\n")

# ============================================================
# PLOT 4 — Wealth fan chart (key strategies only)
# ============================================================
key_strategies <- c("60/40", "MV_plain", "RF", "HybridNaive")
year_full <- 0:N_YEARS

fan_df <- lapply(key_strategies, function(s) {
  wp_mat <- wp[[s]]
  data.frame(
    strategy = s, year = year_full,
    p10 = apply(wp_mat, 2, quantile, 0.10),
    p25 = apply(wp_mat, 2, quantile, 0.25),
    p50 = apply(wp_mat, 2, quantile, 0.50),
    p75 = apply(wp_mat, 2, quantile, 0.75),
    p90 = apply(wp_mat, 2, quantile, 0.90)
  )
}) %>% bind_rows() %>%
  mutate(strategy = factor(strategy, levels = key_strategies),
         family   = strategy_family[as.character(strategy)])

p4 <- ggplot(fan_df, aes(x = year)) +
  geom_ribbon(aes(ymin = p10, ymax = p90, fill = family),
              alpha = 0.25) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = family),
              alpha = 0.40) +
  geom_line(aes(y = p50, color = family), linewidth = 0.7) +
  facet_wrap(~ strategy, ncol = 2, scales = "fixed") +
  scale_fill_manual(values = family_colors) +
  scale_color_manual(values = family_colors) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 42)) +
  labs(
    title    = "Wealth accumulation paths (selected strategies)",
    subtitle = "Wealth normalised to annual salary at age 25. Dark band: 25-75 pct. Light: 10-90 pct.",
    x        = "Career year",
    y        = "Wealth (multiples of initial salary)",
    color    = "Strategy family",
    fill     = "Strategy family",
    caption  = "Source: simulation on DNB CP2022 scenarios."
  ) +
  theme_thesis()

ggsave("Output/Plots/04_wealth_fan.pdf", p4,
       width = 9, height = 7, device = "pdf")
cat("Saved 04_wealth_fan.pdf\n")

# ============================================================
# PLOT 5 — Trade-off scatter
# ============================================================
tradeoff_df <- metrics %>%
  select(strategy, RR_median, P_shortfall_70pct, CVaR_5pct, MDD_mean) %>%
  mutate(strategy = factor(strategy, levels = strategy_order),
         family   = strategy_family[as.character(strategy)])

p5 <- ggplot(tradeoff_df,
             aes(x = P_shortfall_70pct, y = RR_median, color = family)) +
  geom_point(aes(size = -MDD_mean), alpha = 0.85) +
  geom_text(aes(label = strategy), vjust = -1.2, size = 3.2,
            color = "black", show.legend = FALSE) +
  scale_color_manual(values = family_colors) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.65, 1.0)) +
  scale_y_continuous(limits = c(0.42, 0.58),
                     breaks = seq(0.40, 0.60, 0.025),
                     labels = label_number(accuracy = 0.01)) +
  scale_size_continuous(name = "Mean max drawdown",
                        range = c(3, 10),
                        breaks = c(0.12, 0.16, 0.20),
                        labels = c("-12%", "-16%", "-20%")) +
  labs(
    title    = "Risk-return trade-off across the nine strategies",
    subtitle = "Upper-left is preferable: high median RR, low shortfall.",
    x        = "P(RR < 0.70) - Shortfall probability",
    y        = "Median Replacement Ratio",
    color    = "Strategy family",
    caption  = "Point size scales with maximum drawdown (larger = worse drawdown)."
  ) +
  theme_thesis()

ggsave("Output/Plots/05_tradeoff_scatter.pdf", p5,
       width = 9, height = 6, device = "pdf")
cat("Saved 05_tradeoff_scatter.pdf\n")

# ============================================================
# PLOT 6 — Hybrid composition
# ============================================================
comp_df <- data.frame(
  year       = year_seq,
  c_6040     = colMeans(w_6040),
  c_MV_plain = colMeans(w_mv_p),
  c_RF       = colMeans(w_rf),
  c_Hybrid   = colMeans(w_hybrid)
) %>%
  pivot_longer(-year, names_to = "component", values_to = "weight") %>%
  mutate(component = case_when(
    component == "c_6040"     ~ "60/40 component",
    component == "c_MV_plain" ~ "MV plain component",
    component == "c_RF"       ~ "RF component",
    component == "c_Hybrid"   ~ "HybridNaive (combined)"
  )) %>%
  mutate(component = factor(component, levels = c(
    "60/40 component", "MV plain component", "RF component",
    "HybridNaive (combined)"
  )))

comp_colors <- c(
  "60/40 component"        = "#1f77b4",
  "MV plain component"     = "#2ca02c",
  "RF component"           = "#d62728",
  "HybridNaive (combined)" = "#9467bd"
)
comp_linetype <- c(
  "60/40 component"        = "dotted",
  "MV plain component"     = "dotted",
  "RF component"           = "dotted",
  "HybridNaive (combined)" = "solid"
)

p6 <- ggplot(comp_df,
             aes(x = year, y = weight, color = component, linetype = component)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = comp_colors) +
  scale_linetype_manual(values = comp_linetype) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.35, 0.75)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 42)) +
  labs(
    title    = "Naive Hybrid: composition over time",
    subtitle = "Hybrid = (1/3) x 60/40 + (1/3) x MV plain + (1/3) x RF",
    x        = "Career year",
    y        = "Mean equity weight across scenarios",
    color    = NULL, linetype = NULL,
    caption  = paste(
      "60/40 and MV are static. RF varies because it conditions on the realised macro regime.",
      "\nHybrid inherits RF's variation, dampened by the two static components."
    )
  ) +
  theme_thesis()

ggsave("Output/Plots/06_hybrid_composition.pdf", p6,
       width = 8, height = 5, device = "pdf")
cat("Saved 06_hybrid_composition.pdf\n")

cat("\nAll 6 plots saved to Output/Plots/\n")
cat("\nSummary:\n")
cat("  01: Density facet of RR per strategy\n")
cat("  02: Shortfall ranking bar chart\n")
cat("  03: Equity allocation evolution\n")
cat("  04: Wealth fan chart (4 key strategies)\n")
cat("  05: Risk-return trade-off scatter\n")
cat("  06: Hybrid composition (deep dive)\n")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
