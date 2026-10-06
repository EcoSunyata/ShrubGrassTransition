# ==============================================================================
# PNAS Supplementary Information
# Trait classification via k-means + PCA biplots + regression analyses
# Author: Jingyao Sun
# Data source: Dataset_S1.xlsx
# Usage: Delete the first line of the excel sheet before use, since it contains unit info
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(MASS)
library(patchwork)

# ------------------------------------------------------------------------------
# 2. Publication theme
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
data_trait <- read.xlsx("Dataset_S1.xlsx", sheet = "Date_trait")



batch <- data_trait %>%
  filter(vegetation_type %in% c("natural", "boundary")) %>%
  mutate(
    height      = log(height),
    crown       = log(crown),
    base        = log(base),
    SLA         = log(leaf_area / leaf_weight),
    leaf_area   = log(leaf_area),
    leaf_weight = log(leaf_weight),
    leaf_thick  = log(leaf_thick),
    leaf_width  = log(leaf_width),
    leaf_length = log(leaf_length),
    CN          = log(TC / TN),
    TN          = log(TN),
    TC          = log(TC),
    TP          = log(TP)
  ) %>%
  drop_na(height, crown, base, SLA,
          leaf_area, leaf_weight, leaf_thick, leaf_width, leaf_length,
          TN, TC, TP)

# ------------------------------------------------------------------------------
# 4. PCA on functional traits
# ------------------------------------------------------------------------------
pca_result <- prcomp(
  ~ height + crown + base + SLA +
    leaf_area + leaf_weight + leaf_thick + leaf_width + leaf_length +
    TN + TC + TP + CN,
  data = batch,
  center = TRUE,
  scale = TRUE
)

fit <- summary(pca_result)
variance_table <- fit$importance %>% t() %>% as.data.frame()

loadings <- pca_result$rotation[, 1:10] * 6
loadings <- as.data.frame(loadings)
loadings$variable <- rownames(loadings)

# ------------------------------------------------------------------------------
# 5. Project traits into PCA space
# ------------------------------------------------------------------------------
pca_data <- pca_result$x %>% as.data.frame()
pca_data$precipitation   <- batch$precipitation
pca_data$class_species   <- batch$class_species
pca_data$id              <- batch$id
pca_data$cover           <- batch$cover
pca_data$totalcover      <- batch$totalcover
pca_data$code            <- batch$code
pca_data$vegetation_type <- batch$vegetation_type
pca_data$species         <- batch$species

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
# 7. Dominant species class per plot
# ------------------------------------------------------------------------------
dominant_class <- pca_data %>%
  group_by(code, class_species) %>%
  summarise(total_cover = sum(cover, na.rm = TRUE), .groups = "drop") %>%
  group_by(code) %>%
  slice_max(order_by = total_cover, n = 1) %>%
  dplyr::select(code, dominant_class = class_species)

batch_summary <- batch_summary %>%
  left_join(dominant_class, by = "code")

# ------------------------------------------------------------------------------
# 8. PCA axis label helper
# ------------------------------------------------------------------------------
get_pca_label <- function(axis) {
  paste0(axis, " (",
         round(variance_table[axis, "Proportion of Variance"] * 100, 1),
         "%)")
}

# ==============================================================================
# 9. Automatic trait classification via k-means (k = 2)
# ==============================================================================
loadings_matrix <- as.matrix(loadings[, c("PC2", "PC3")])
rownames(loadings_matrix) <- loadings$variable

set.seed(123)
km <- kmeans(loadings_matrix, centers = 2, nstart = 25)
loadings$cluster <- km$cluster

# Anchor traits to assign biological meaning
anchor_acq  <- "SLA"
anchor_cons <- "leaf_weight"

cluster_acq  <- loadings$cluster[loadings$variable == anchor_acq]
cluster_cons <- loadings$cluster[loadings$variable == anchor_cons]

if (cluster_acq == cluster_cons) {
  warning("Anchor traits fell into the same cluster — clustering failed to separate them.")
}

loadings$auto_class <- ifelse(loadings$cluster == cluster_acq,
                              "Acquisitive", "Conservative")

conservative_vars <- loadings$variable[loadings$auto_class == "Conservative"]
acquisitive_vars  <- loadings$variable[loadings$auto_class == "Acquisitive"]

cat("Conservative traits:", paste(conservative_vars, collapse = ", "), "\n")
cat("Acquisitive traits: ", paste(acquisitive_vars,  collapse = ", "), "\n")

# ==============================================================================
# 10. Figure: Automatic trait classification
# ==============================================================================
p_cluster <- ggplot(loadings, aes(x = PC2, y = PC3)) +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(aes(color = auto_class), size = 4.5, alpha = 0.9) +
  geom_text(aes(label = variable), vjust = -1.2, size = 3.2, color = "black") +
  scale_color_manual(values = c("Conservative" = "#2C7FB8",
                                "Acquisitive"  = "#D95F02"),
                     name = "Trait class") +
  labs(x = paste0("PC2 (",
                  round(variance_table["PC2", "Proportion of Variance"] * 100, 1), "%)"),
       y = paste0("PC3 (",
                  round(variance_table["PC3", "Proportion of Variance"] * 100, 1), "%)"),
       title = "Automatic trait classification (k-means, k = 2)") +
  theme_bw() + mytheme

print(p_cluster)
#ggsave("Fig_trait_clustering.pdf", p_cluster, width = 5.5, height = 4.5)

# ==============================================================================
# 11. Define continuous trait axes (PC1–PC2 and PC2–PC3)
# ==============================================================================
v_conservative_123 <- c(
  mean(loadings$PC1[loadings$auto_class == "Conservative"]),
  mean(loadings$PC2[loadings$auto_class == "Conservative"]),
  mean(loadings$PC3[loadings$auto_class == "Conservative"])
)
v_acquisitive_123 <- c(
  mean(loadings$PC1[loadings$auto_class == "Acquisitive"]),
  mean(loadings$PC2[loadings$auto_class == "Acquisitive"]),
  mean(loadings$PC3[loadings$auto_class == "Acquisitive"])
)

v_trait_123 <- v_acquisitive_123 - v_conservative_123
v_trait_123 <- v_trait_123 / sqrt(sum(v_trait_123^2))

cat("\nContinuous trait axis (PC1, PC2, PC3):", round(v_trait_123, 3), "\n")

v_trait_12 <- v_trait_123[c(1, 2)]
v_trait_12 <- v_trait_12 / sqrt(sum(v_trait_12^2))

v_trait_23 <- v_trait_123[c(2, 3)]
v_trait_23 <- v_trait_23 / sqrt(sum(v_trait_23^2))

# ==============================================================================
# 12. Prepare plot-level data
# ==============================================================================
dat <- batch_summary %>%
  filter(vegetation_type == "natural") %>%
  dplyr::select(code, PC1, PC2, PC3, size_mean, distance_mean, dominant_class) %>%
  drop_na(size_mean, distance_mean)

# ==============================================================================
# 13. Patch-characteristic direction
#     Combines patch size AND nearest-neighbor distance into a single axis
# ==============================================================================

# --- 13.1 Combine patch size and distance into a single composite ---
# PCA on the two patch metrics to obtain a single "patch characteristic" axis
patch_pca <- prcomp(~ size_mean + distance_mean,
                    data = dat, center = TRUE, scale. = TRUE)

dat$patch_composite <- patch_pca$x[, 1]

# Ensure the composite is oriented positively:
# higher values = larger patches and larger spacing
if (cor(dat$patch_composite, dat$size_mean, use = "complete.obs") < 0) {
  dat$patch_composite <- -dat$patch_composite
}

# Report loadings of the patch PCA
cat("\nPatch PCA loadings:\n")
print(round(patch_pca$rotation, 3))

# --- 13.2 Regress PC scores on the composite patch characteristic ---
fit_PC1 <- lm(PC1 ~ patch_composite, data = dat)
fit_PC2 <- lm(PC2 ~ patch_composite, data = dat)
fit_PC3 <- lm(PC3 ~ patch_composite, data = dat)

v_patch_123 <- c(coef(fit_PC1)[2], coef(fit_PC2)[2], coef(fit_PC3)[2])
v_patch_123 <- v_patch_123 / sqrt(sum(v_patch_123^2))

cat("\nPatch-characteristic direction (PC1, PC2, PC3):",
    round(v_patch_123, 3), "\n")

# 2D projections
v_patch_12 <- v_patch_123[c(1, 2)]
v_patch_12 <- v_patch_12 / sqrt(sum(v_patch_12^2))

v_patch_23 <- v_patch_123[c(2, 3)]
v_patch_23 <- v_patch_23 / sqrt(sum(v_patch_23^2))

# ==============================================================================
# 14. Project plots onto each plane's trait axis
# ==============================================================================
dat <- dat %>%
  mutate(
    proj_trait_12 = PC1 * v_trait_12[1] + PC2 * v_trait_12[2],
    proj_trait_23 = PC2 * v_trait_23[1] + PC3 * v_trait_23[2]
  )

# ==============================================================================
# 15. Regressions in each plane (patch size and distance separately)
# ==============================================================================

# --- PC1–PC2 plane ---
fit_size_12 <- lm(size_mean ~ proj_trait_12, data = dat)
r2_size_12  <- summary(fit_size_12)$r.squared
p_size_12   <- summary(fit_size_12)$coefficients["proj_trait_12", "Pr(>|t|)"]

fit_dist_12 <- lm(distance_mean ~ proj_trait_12, data = dat)
r2_dist_12  <- summary(fit_dist_12)$r.squared
p_dist_12   <- summary(fit_dist_12)$coefficients["proj_trait_12", "Pr(>|t|)"]

# --- PC2–PC3 plane ---
fit_size_23 <- lm(size_mean ~ proj_trait_23, data = dat)
r2_size_23  <- summary(fit_size_23)$r.squared
p_size_23   <- summary(fit_size_23)$coefficients["proj_trait_23", "Pr(>|t|)"]

fit_dist_23 <- lm(distance_mean ~ proj_trait_23, data = dat)
r2_dist_23  <- summary(fit_dist_23)$r.squared
p_dist_23   <- summary(fit_dist_23)$coefficients["proj_trait_23", "Pr(>|t|)"]

cat("\n--- PC1–PC2 plane ---\n")
cat("Size    : R² =", round(r2_size_12, 4),
    " p =", format(p_size_12, scientific = TRUE), "\n")
cat("Distance: R² =", round(r2_dist_12, 4),
    " p =", format(p_dist_12, scientific = TRUE), "\n")

cat("\n--- PC2–PC3 plane ---\n")
cat("Size    : R² =", round(r2_size_23, 4),
    " p =", format(p_size_23, scientific = TRUE), "\n")
cat("Distance: R² =", round(r2_dist_23, 4),
    " p =", format(p_dist_23, scientific = TRUE), "\n")

# ==============================================================================
# 16. Figure A: PC1–PC2 PCA biplot with arrows
# ==============================================================================
arrow_scale <- 3

p_pc1_pc2 <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(data = pca_data %>% filter(class_species == "shrub"),
             aes(x = PC1, y = PC2), color = "#8DA0CB", alpha = 0.2, size = 1.2) +
  geom_point(data = pca_data %>% filter(class_species == "herb"),
             aes(x = PC1, y = PC2), color = "#FC8D62", alpha = 0.2, size = 1.2) +
  geom_point(data = dat, aes(x = PC1, y = PC2, color = size_mean),
             size = 3.5, shape = 15, alpha = 0.9) +
  annotate("segment", x = 0, y = 0,
           xend = v_trait_12[1] * arrow_scale,
           yend = v_trait_12[2] * arrow_scale,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#2C7FB8", size = 1.4) +
  annotate("text",
           x = v_trait_12[1] * arrow_scale * 1.2,
           y = v_trait_12[2] * arrow_scale * 1.2,
           label = "Trait axis\n(conservative → acquisitive)",
           color = "#2C7FB8", size = 3, fontface = "bold", hjust = 0.5) +
  annotate("segment", x = 0, y = 0,
           xend = v_patch_12[1] * arrow_scale,
           yend = v_patch_12[2] * arrow_scale,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#D95F02", size = 1.4) +
  annotate("text",
           x = v_patch_12[1] * arrow_scale * 1.2,
           y = v_patch_12[2] * arrow_scale * 1.2,
           label = "Patch-characteristic\ndirection\n(size + distance)",
           color = "#D95F02", size = 3, fontface = "bold", hjust = 0.5) +
  scale_color_gradient2(low = "#f8ef20", mid = "#1f918d", high = "#44035b",
                        midpoint = median(dat$size_mean, na.rm = TRUE),
                        name = "Patch size") +
  labs(x = get_pca_label("PC1"), y = get_pca_label("PC2")) +
  theme_bw() + mytheme +
  theme(panel.grid = element_blank(), legend.position = "right")

print(p_pc1_pc2)
#ggsave("Fig_PC1_PC2.pdf", p_pc1_pc2, width = 5.5, height = 4.5)

# ==============================================================================
# 17. Figure B: PC2–PC3 PCA biplot with arrows
# ==============================================================================
p_pc2_pc3 <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey") +
  geom_point(data = pca_data %>% filter(class_species == "shrub"),
             aes(x = PC2, y = PC3), color = "#8DA0CB", alpha = 0.2, size = 1.2) +
  geom_point(data = pca_data %>% filter(class_species == "herb"),
             aes(x = PC2, y = PC3), color = "#FC8D62", alpha = 0.2, size = 1.2) +
  geom_point(data = dat, aes(x = PC2, y = PC3, color = size_mean),
             size = 3.5, shape = 15, alpha = 0.9) +
  annotate("segment", x = 0, y = 0,
           xend = v_trait_23[1] * arrow_scale,
           yend = v_trait_23[2] * arrow_scale,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#2C7FB8", size = 1.4) +
  annotate("text",
           x = v_trait_23[1] * arrow_scale * 1.2,
           y = v_trait_23[2] * arrow_scale * 1.2,
           label = "Trait axis\n(conservative → acquisitive)",
           color = "#2C7FB8", size = 3, fontface = "bold", hjust = 0.5) +
  annotate("segment", x = 0, y = 0,
           xend = v_patch_23[1] * arrow_scale,
           yend = v_patch_23[2] * arrow_scale,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#D95F02", size = 1.4) +
  annotate("text",
           x = v_patch_23[1] * arrow_scale * 1.2,
           y = v_patch_23[2] * arrow_scale * 1.2,
           label = "Patch-characteristic\ndirection\n(size + distance)",
           color = "#D95F02", size = 3, fontface = "bold", hjust = 0.5) +
  scale_color_gradient2(low = "#f8ef20", mid = "#1f918d", high = "#44035b",
                        midpoint = median(dat$size_mean, na.rm = TRUE),
                        name = "Patch size") +
  labs(x = get_pca_label("PC2"), y = get_pca_label("PC3")) +
  theme_bw() + mytheme +
  theme(panel.grid = element_blank(), legend.position = "right")

print(p_pc2_pc3)
#ggsave("Fig_PC2_PC3.pdf", p_pc2_pc3, width = 5.5, height = 4.5)

# ==============================================================================
# 18. Figure C: Regression in PC1–PC2 plane
# ==============================================================================
dat_long_12 <- dat %>%
  mutate(
    size_z     = scale(size_mean)[, 1],
    distance_z = scale(distance_mean)[, 1]
  ) %>%
  pivot_longer(cols = c(size_z, distance_z),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
                           "size_z"     = "Patch size",
                           "distance_z" = "Nearest-neighbor distance"))

p_reg_12 <- ggplot(dat_long_12,
                   aes(x = proj_trait_12, y = value, color = variable)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
  scale_color_manual(values = c("Patch size" = "#2C7FB8",
                                "Nearest-neighbor distance" = "#D95F02"),
                     name = NULL) +
  annotate("text", x = min(dat$proj_trait_12), y = max(dat_long_12$value),
           label = sprintf("Patch size: R² = %.3f, p = %.2e",
                           r2_size_12, p_size_12),
           hjust = 0, vjust = 1, size = 3, color = "#2C7FB8") +
  annotate("text", x = min(dat$proj_trait_12), y = max(dat_long_12$value) - 0.4,
           label = sprintf("Distance: R² = %.3f, p = %.2e",
                           r2_dist_12, p_dist_12),
           hjust = 0, vjust = 1, size = 3, color = "#D95F02") +
  labs(x = "Trait score in PC1–PC2 plane (conservative → acquisitive)",
       y = "Standardized patch characteristic (z-score)") +
  theme_bw() + mytheme

print(p_reg_12)
#ggsave("Fig_reg_PC1_PC2.pdf", p_reg_12, width = 6, height = 4.5)

# ==============================================================================
# 19. Figure D: Regression in PC2–PC3 plane
# ==============================================================================
dat_long_23 <- dat %>%
  mutate(
    size_z     = scale(size_mean)[, 1],
    distance_z = scale(distance_mean)[, 1]
  ) %>%
  pivot_longer(cols = c(size_z, distance_z),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
                           "size_z"     = "Patch size",
                           "distance_z" = "Nearest-neighbor distance"))

p_reg_23 <- ggplot(dat_long_23,
                   aes(x = proj_trait_23, y = value, color = variable)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
  scale_color_manual(values = c("Patch size" = "#2C7FB8",
                                "Nearest-neighbor distance" = "#D95F02"),
                     name = NULL) +
  annotate("text", x = min(dat$proj_trait_23), y = max(dat_long_23$value),
           label = sprintf("Patch size: R² = %.3f, p = %.2e",
                           r2_size_23, p_size_23),
           hjust = 0, vjust = 1, size = 3, color = "#2C7FB8") +
  annotate("text", x = min(dat$proj_trait_23), y = max(dat_long_23$value) - 0.4,
           label = sprintf("Distance: R² = %.3f, p = %.2e",
                           r2_dist_23, p_dist_23),
           hjust = 0, vjust = 1, size = 3, color = "#D95F02") +
  labs(x = "Trait score in PC2–PC3 plane (conservative → acquisitive)",
       y = "Standardized patch characteristic (z-score)") +
  theme_bw() + mytheme

print(p_reg_23)
#ggsave("Fig_reg_PC2_PC3.pdf", p_reg_23, width = 6, height = 4.5)

