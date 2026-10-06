# ==============================================================================
# PNAS Supplementary Information
# Partial effect analysis of soil and topographic variables on vegetation cover
# Author: Jingyao Sun
# Data source: Dataset_S1.xlsx
# Usage: Delete the first line of the excel sheet before use, since it contains unit info
# ==============================================================================

library(tidyverse)
library(broom)
library(openxlsx)
library(ggplot2)

data_analysis <- read.xlsx("Dataset_S1.xlsx", sheet = "Data_site")


# ============================================================
# 1. Partial effect analysis: independent effects of soil factors
#    after controlling for precipitation
# ============================================================
soil_vars <- c(
  "elevation",
  "BD_0-5cm", "CEC_0-5cm", "pH_0-5cm", "clay_0-5cm", "silt_0-5cm", "sand_0-5cm",
  "OC_0-5cm", "TN_0-5cm", "TP_0-5cm",
  "BD_5-15cm", "CEC_5-15cm", "pH_5-15cm", "clay_5-15cm", "silt_5-15cm", "sand_5-15cm",
  "OC_5-15cm", "TN_5-15cm", "TP_5-15cm",
  "BD_15-30cm", "CEC_15-30cm", "pH_15-30cm", "clay_15-30cm", "silt_15-30cm", "sand_15-30cm",
  "OC_15-30cm", "TN_15-30cm", "TP_15-30cm",
  "BD_30-60cm", "CEC_30-60cm", "pH_30-60cm", "clay_30-60cm", "silt_30-60cm", "sand_30-60cm",
  "OC_30-60cm", "TN_30-60cm", "TP_30-60cm",
  "BD_60-100cm", "CEC_60-100cm", "pH_60-100cm", "clay_60-100cm", "silt_60-100cm", "sand_60-100cm",
  "OC_60-100cm", "TN_60-100cm", "TP_60-100cm"
)

# For each soil factor, fit cover ~ precipitation + soil_var
# Extract the partial effect coefficient, standard error, 95% CI, and p-value
partial_effects <- lapply(soil_vars, function(var) {
  form <- as.formula(paste0("shrub ~ precipitation + `", var, "`"))  # cover/shrub/grass
  mod <- lm(form, data = data_analysis)
  
  coef_tab <- tidy(mod, conf.int = TRUE) %>%
    filter(term == paste0("`", var, "`"))
  
  coef_tab$variable <- var
  coef_tab
}) %>% bind_rows()

# Output results
partial_effects %>%
  dplyr::select(variable, estimate, conf.low, conf.high, p.value) %>%
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  arrange(p.value) %>%
  print(n = 30) %>% view()

# ============================================================
# 3. Forest plot: partial effect coefficients + 95% CI
# ============================================================
if (!"elevation" %in% partial_effects$variable) {
  fit_elev <- lm(cover ~ precipitation + elevation, data = data_analysis)
  elev_effect <- tidy(fit_elev, conf.int = TRUE) %>%
    filter(term == "elevation") %>%
    mutate(variable = "elevation")
  partial_effects <- bind_rows(partial_effects, elev_effect)
}

partial_effects_plot <- partial_effects %>%
  mutate(
    is_elevation = variable == "elevation",
    indicator = ifelse(is_elevation, "elevation", sub("_.*$", "", variable)),
    depth = ifelse(is_elevation, "", sub("^[^_]*_", "", variable)),
    
    p_label = case_when(
      p.value < 0.001 ~ paste0("p=", format(p.value, scientific = TRUE, digits = 2)),
      p.value < 0.05  ~ paste0("p=", sprintf("%.3f", p.value)),
      TRUE            ~ NA_character_
    ),
    sig_class = ifelse(p.value < 0.05, "Significant", "Not significant"),
    
    indicator_order = match(indicator,
                            c("sand", "silt", "clay",
                              "BD", "CEC", "pH",
                              "OC", "TN", "TP",
                              "elevation")),
    depth_order = match(depth, c("0-5cm", "5-15cm", "15-30cm",
                                 "30-60cm", "60-100cm", "")),
    label = ifelse(is_elevation,
                   "elevation",
                   paste0(indicator, " (", depth, ")")),
    sort_key = indicator_order * 100 + depth_order
  ) %>%
  arrange(sort_key) %>%
  mutate(label = factor(label, levels = unique(label)))

# ============================================================
# 3. Label x position: close to each conf.high
# ============================================================
partial_effects_plot <- partial_effects_plot %>%
  mutate(
    offset = 0.02 * diff(range(c(conf.low, conf.high), na.rm = TRUE)),
    label_x = conf.high + offset
  )

# ============================================================
# 4. Plot
# ============================================================
p <- ggplot(partial_effects_plot, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, color = "grey30") +
  geom_point(aes(color = sig_class), size = 2.5) +
  geom_text(aes(x = label_x, label = p_label),
            hjust = 0, size = 2.4, color = "black") +
  scale_color_manual(values = c("Significant" = "#D95F02",
                                "Not significant" = "#2C7FB8"),
                     name = NULL) +
  labs(x = "Partial effect on cover (controlling for precipitation)",
       y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey95"),
    legend.position = "bottom",
    axis.text.y = element_text(size = 7)
  ) +
  # Only extend the right side a little, not too much
  coord_cartesian(
    xlim = c(min(partial_effects_plot$conf.low, na.rm = TRUE) * 1.1,
             max(partial_effects_plot$label_x, na.rm = TRUE) * 1.05)
  )

print(p)

#ggsave("Fig_partial_effects.pdf", p, width = 3, height = 6)