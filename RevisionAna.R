# ==============================================================================
# PNAS Supplementary Information
# Analysis of plant functional traits, patch characteristics, and environmental covariates
# 
# Author: Jingyao Sun
# Description: This script performs PCA on plant functional traits, maps patch 
#   characteristics onto trait space, and tests the independent contribution of 
#   environmental covariates to vegetation cover.
# 
# Data source: DataAppendixAna.xlsx
# Output: PCA ordination plots, trait–patch relationship figures, and regression summaries
# Ussage: Delete the first line of the excel sheet before use, since it contains unit info
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(broom)
library(openxlsx)

# ------------------------------------------------------------------------------
# 2. Define publication theme
# ------------------------------------------------------------------------------

mytheme <- theme(
  panel.background = element_rect(fill = "white", color = NA),
  panel.border = element_rect(color = "black", fill = NA, size = 0.8),
  panel.grid = element_blank(),
  axis.text = element_text(color = "black", size = 8),
  axis.title = element_text(color = "black", size = 10),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  legend.key.size = unit(0.3, "cm"),
  strip.background = element_blank(),
  strip.text = element_text(size = 9)
)

# ------------------------------------------------------------------------------
# 3. Load and transform trait data
# ------------------------------------------------------------------------------

data_trait <- read.xlsx("Dataset_S1.xlsx", sheet = "Data_trait")

batch <- data_trait %>%
  filter(vegetation_type %in% c("natural", "boundary")) %>%
  mutate(
    height = log(height),
    crown = log(crown),
    base = log(base),
    SLA = log(leaf_area / leaf_weight),
    leaf_area = log(leaf_area),
    leaf_weight = log(leaf_weight),
    leaf_thick = log(leaf_thick),
    leaf_width = log(leaf_width),
    leaf_length = log(leaf_length),
    CN = log(TC / TN),
    TN = log(TN),
    TC = log(TC),
    TP = log(TP)
  ) %>%
  drop_na(height, crown, base, SLA,
          leaf_area, leaf_weight, leaf_thick, leaf_width, leaf_length,
          TN, TC, TP)

# ------------------------------------------------------------------------------
# 4. PCA
# ------------------------------------------------------------------------------

pca_result <- prcomp(
  ~ height + crown + base + SLA +
    leaf_area + leaf_weight + leaf_thick + leaf_width + leaf_length +
    TN + TC + TP + CN,
  data = batch,
  center = TRUE,
  scale = TRUE
)

# Variance explained
fit <- summary(pca_result)
variance_table <- fit$importance %>% t() %>% as.data.frame()

# Loadings (scaled for visualization)
loadings <- pca_result$rotation[, 1:10] * 6
loadings <- as.data.frame(loadings)
loadings$variable <- rownames(loadings)

# ------------------------------------------------------------------------------
# 5. Project traits into PCA space
# ------------------------------------------------------------------------------

pca_data <- pca_result$x %>% as.data.frame()
pca_data$precipitation <- batch$precipitation
pca_data$class_species <- batch$class_species
pca_data$id <- batch$id
pca_data$cover <- batch$cover
pca_data$totalcover <- batch$totalcover
pca_data$code <- batch$code
pca_data$vegetation_type <- batch$vegetation_type
pca_data$species <- batch$species

# ------------------------------------------------------------------------------
# 6. Cover-weighted mean trait scores per plot
# ------------------------------------------------------------------------------

data_patch <- read.xlsx("Dataset_S1.xlsx", sheet = "Data_site")

pca_data <- pca_data %>%
  mutate(
    weight = cover / totalcover,
    PC1w = PC1 * weight,
    PC2w = PC2 * weight,
    PC3w = PC3 * weight,
    PC4w = PC4 * weight
  )

trait_summary <- pca_data %>%
  group_by(code, vegetation_type) %>%
  summarise(
    PC1 = sum(PC1w),
    PC2 = sum(PC2w),
    PC3 = sum(PC3w),
    PC4 = sum(PC4w),
    .groups = "drop"
  )

batch_summary <- merge(
  trait_summary,
  data_patch,
  by = c("code", "vegetation_type")
)

# ------------------------------------------------------------------------------
# 7. PCA biplots
# ------------------------------------------------------------------------------

# Helper for axis labels
get_pca_label <- function(axis) {
  paste0(axis, " (", round(variance_table[axis, "Proportion of Variance"] * 100, 1), "%)")
}

# ---- Fig. 2I: PC1 vs PC2, colored by patch size ----
p_pc1_pc2_size <- batch_summary %>%
  filter(vegetation_type == "natural") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(data = pca_data %>% filter(class_species == "shrub"),
             aes(x = PC1, y = PC2), color = "#8DA0CB", alpha = 0.6) +
  geom_point(data = pca_data %>% filter(class_species == "herb"),
             aes(x = PC1, y = PC2), color = "#FC8D62", alpha = 0.6) +
  geom_point(aes(color = size_mean), size = 3, shape = 15) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgrey") +
  geom_text(data = loadings, aes(label = variable), size = 3) +
  scale_color_gradient2(low = "#f8ef20", mid = "#1f918d", high = "#44035b",
                        midpoint = 1, name = "Patch size") +
  labs(x = get_pca_label("PC1"), y = get_pca_label("PC2")) +
  theme_bw() + mytheme +
  theme(panel.grid = element_blank(), legend.position = "right")
print(p_pc1_pc2_size)
# ---- Fig. 2J: PC2 vs PC3, colored by patch size ----
p_pc2_pc3_size <- batch_summary %>%
  filter(vegetation_type == "natural") %>%
  ggplot(aes(x = PC2, y = PC3)) +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(data = pca_data %>% filter(class_species == "shrub"),
             aes(x = PC2, y = PC3), color = "#8DA0CB", alpha = 0.6) +
  geom_point(data = pca_data %>% filter(class_species == "herb"),
             aes(x = PC2, y = PC3), color = "#FC8D62", alpha = 0.6) +
  geom_point(aes(color = size_mean), size = 3, shape = 15) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgrey") +
  geom_text(data = loadings, aes(label = variable), size = 3) +
  scale_color_gradient2(low = "#f8ef20", mid = "#1f918d", high = "#44035b",
                        midpoint = 0.8, name = "Patch size") +
  labs(x = get_pca_label("PC2"), y = get_pca_label("PC3")) +
  theme_bw() + mytheme +
  theme(panel.grid = element_blank(), legend.position = "right")
print(p_pc2_pc3_size)
# ---- PC1 vs PC2, colored by nearest-neighbor distance ----
p_pc1_pc2_dist <- batch_summary %>%
  filter(vegetation_type == "natural") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(data = pca_data %>% filter(class_species == "shrub"),
             aes(x = PC1, y = PC2), color = "#8DA0CB", alpha = 0.6) +
  geom_point(data = pca_data %>% filter(class_species == "herb"),
             aes(x = PC1, y = PC2), color = "#FC8D62", alpha = 0.6) +
  geom_point(aes(color = distance_mean), size = 3, shape = 15) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgrey") +
  geom_text(data = loadings, aes(label = variable), size = 3) +
  scale_color_gradient2(low = "#f8ef20", mid = "#1f918d", high = "#44035b",
                        midpoint = 2, name = "Patch distance") +
  labs(x = get_pca_label("PC1"), y = get_pca_label("PC2")) +
  theme_bw() + mytheme +
  theme(panel.grid = element_blank(), legend.position = "right")
print(p_pc1_pc2_dist)
# ---- PC2 vs PC3, colored by nearest-neighbor distance ----
p_pc2_pc3_dist <- batch_summary %>%
  filter(vegetation_type == "natural") %>%
  ggplot(aes(x = PC2, y = PC3)) +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(data = pca_data %>% filter(class_species == "shrub"),
             aes(x = PC2, y = PC3), color = "#8DA0CB", alpha = 0.6) +
  geom_point(data = pca_data %>% filter(class_species == "herb"),
             aes(x = PC2, y = PC3), color = "#FC8D62", alpha = 0.6) +
  geom_point(aes(color = distance_mean), size = 3, shape = 15) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgrey") +
  geom_text(data = loadings, aes(label = variable), size = 3) +
  scale_color_gradient2(low = "#f8ef20", mid = "#1f918d", high = "#44035b",
                        midpoint = 2, name = "Patch distance") +
  labs(x = get_pca_label("PC2"), y = get_pca_label("PC3")) +
  theme_bw() + mytheme +
  theme(panel.grid = element_blank(), legend.position = "right")
print(p_pc2_pc3_dist)
# Save if needed
# ggsave("Fig2I.pdf", p_pc1_pc2_size, width = 6, height = 5)
# ggsave("Fig2J.pdf", p_pc2_pc3_size, width = 6, height = 5)

# ------------------------------------------------------------------------------
# 8. Covariate analysis: cover ~ precipitation + soil_factor
# ------------------------------------------------------------------------------

data_analysis <- read.xlsx("Dataset_S1.xlsx", sheet = "Data_site")

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

# Run all two-predictor regressions
regression_results <- lapply(soil_vars, function(var) {
  # could change "cover" to "shrub"/"grass"
  form <- as.formula(paste0("cover ~ precipitation + `", var, "`")) 
  mod <- lm(form, data = data_analysis)
  tidy(mod) %>% mutate(variable = var)
}) %>% bind_rows()

# Show significant predictors only
regression_results %>%
  filter(p.value < 0.05 & term != "(Intercept)") %>%
  view()

