##########################################################
# Proliferative Module Score Analysis in T Cells
# Author: Qing Zhong
##########################################################

library(tidyverse)
library(Seurat)
library(ggpubr)

# === 1. Load Data ===
load("step04_Major_Annotated.RData")
T_Cell <- readRDS("step05_Annotated_Tcell_Integrated.RDS")

# === 2. Identify Proliferative Signature Genes ===
Proliferative_Signature <- FindMarkers(
  T_Cell,
  ident.1 = c("CD8_06_Proliferative_MKI67", "CD4_04_Proliferative_MKI67"),
  group.by = "minor_Cluster"
) %>%
  filter(avg_log2FC > 2 & pct.1 > 0.3) %>%
  rownames()

# === 3. Add Module Score to Seurat Object ===
T_Cell <- AddModuleScore(
  T_Cell,
  features = list(Proliferative_Signature),
  name = "Proliferative"
)

# === 4. Define Function for Module Score Plotting ===
plot_module_score <- function(meta_data, module_col, title, comparison_list, y_limits = NULL, label_y = NULL) {
  meta_data %>%
    mutate(Group = factor(Group, levels = c("Control", "Olaparib", "Concurrent", "AZD1775", "Sequential"))) %>%
    filter(minor_Cluster != "NK") %>%
    mutate(score = .data[[module_col]] - min(.data[[module_col]])) %>%
    ggplot(aes(x = Group, y = score, fill = Group)) +
    stat_summary(fun = mean, geom = "bar") +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.25) +
    scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6")) +
    theme_classic() +
    labs(title = title) +
    ggpubr::stat_compare_means(comparisons = comparison_list, label.y = label_y) +
    {if (!is.null(y_limits)) coord_cartesian(ylim = y_limits) else NULL}
}

# === 5. Define Comparison Groups ===
comparisons <- list(
  c("Control", "Olaparib"),
  c("Control", "Concurrent"),
  c("Control", "AZD1775"),
  c("Control", "Sequential"),
  c("Sequential", "Concurrent")
)

# === 6. Plot Proliferative Module Score (All T cells except NK) ===
p1 <- plot_module_score(
  T_Cell@meta.data,
  module_col = "Proliferative1",
  title = "Proliferative Status",
  comparison_list = comparisons,
  y_limits = c(0, 0.5),
  label_y = c(0.3, 0.33, 0.36, 0.39, 0.28)
)

# === 7. ANOVA and Plot for CD8 T cells ===
cd8_data <- T_Cell@meta.data %>%
  filter(str_detect(minor_Cluster, "CD8")) %>%
  mutate(
    Proliferative1 = Proliferative1 - min(Proliferative1),
    Group = factor(Group, levels = c("Control", "Olaparib", "Concurrent", "AZD1775", "Sequential"))
  )

p_cd8 <- ggplot(cd8_data, aes(x = Group, y = Proliferative1, fill = Group)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.25) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6")) +
  theme_classic() +
  labs(title = "Proliferative Status in CD8 T Cells") +
  ggpubr::stat_compare_means(method = "anova", label.y = 0.45, label = "p.format") +
  coord_cartesian(ylim = c(0, 0.5))

ggsave(filename = "T_Cell_CD8_Proliferative_ANOVA.pdf", plot = p_cd8, width = 8.74 * 3, height = 5.18 * 3, units = "cm")

# ANOVA statistics for CD8 T cells
anova_result_cd8 <- aov(Proliferative1 ~ Group, data = cd8_data)
print(summary(anova_result_cd8))

tukey_result_cd8 <- TukeyHSD(anova_result_cd8)
print(tukey_result_cd8)

# === 8. ANOVA and Plot for all T cells (excluding DP clusters) ===
T_data <- T_Cell@meta.data %>%
  filter(!str_detect(minor_Cluster, "^DP")) %>%
  mutate(
    Proliferative1 = Proliferative1 - min(Proliferative1),
    Group = factor(Group, levels = c("Control", "Olaparib", "Concurrent", "AZD1775", "Sequential"))
  )

p_T <- ggplot(T_data, aes(x = Group, y = Proliferative1, fill = Group)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.25) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6")) +
  theme_classic() +
  labs(title = "Proliferative Status in T Cells") +
  ggpubr::stat_compare_means(method = "anova", label.y = 0.45, label = "p.format") +
  coord_cartesian(ylim = c(0, 0.5))

ggsave(filename = "T_Cell_Proliferative_ANOVA.pdf", plot = p_T, width = 8.74 * 3, height = 5.18 * 3, units = "cm")

# ANOVA statistics for all T cells
anova_result_T <- aov(Proliferative1 ~ Group, data = T_data)
print(summary(anova_result_T))

tukey_result_T <- TukeyHSD(anova_result_T)
print(tukey_result_T)
