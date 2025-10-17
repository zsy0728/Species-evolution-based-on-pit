# Songyao 2025.04.17 - Breakpoints Difference 小提琴图 + 分三组
rm(list = ls())

library(ggsci)
library(ggplot2)
library(R.matlab)
library(openxlsx)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
# Set working directory
outdir = "Y:/sda/results/Evolution_cortical_shape/figure/"

# ------------------ Load Data ------------------
data <- read.xlsx('Y:/sda/results/Evolution_cortical_shape/statistic/saved/breakpoints_and_differentaiation_time.xlsx')
breakpoints <- data[[2]]
species_names <- data[[1]]

DIST <- readMat('Y:/sda/results/Evolution_cortical_shape/statistic/kinship/Species_kinship_DIST.mat')
DIST <- DIST[["DIST"]]

Info <- read.xlsx("Y:/sda/results/Evolution_cortical_shape/data_info/Species_info.xlsx")
primata_list <- Info$Order == "Primata"
names(primata_list) <- Info$Species

# ------------------ Compute pairwise differences ------------------
num_species <- length(breakpoints)
breakpoints_diff <- c()
dist_diff <- c()
species1_vec <- c()
species2_vec <- c()
pair_type <- c()

for (i in 1:(num_species - 1)) {
  for (j in (i + 1):num_species) {
    species1 <- species_names[i]
    species2 <- species_names[j]
    
    is_p1 <- ifelse(species1 %in% names(primata_list), primata_list[[species1]], FALSE)
    is_p2 <- ifelse(species2 %in% names(primata_list), primata_list[[species2]], FALSE)
    
    type <- ifelse(is_p1 & is_p2, "P-P",
                   ifelse(!is_p1 & !is_p2, "N-N", "N-P"))
    
    breakpoints_diff <- c(breakpoints_diff, abs(breakpoints[i] - breakpoints[j]))
    dist_diff <- c(dist_diff, DIST[i, j])
    species1_vec <- c(species1_vec, species1)
    species2_vec <- c(species2_vec, species2)
    pair_type <- c(pair_type, type)
  }
}
# ------------------ Create Data Frame ------------------
plot_data <- data.frame(
  Species1 = species1_vec,
  Species2 = species2_vec,
  Diff = breakpoints_diff,
  Distance = dist_diff,
  PairType = factor(pair_type, levels = c("P-P", "N-P", "N-N"))
)
# 确保分组顺序
plot_data$PairType <- factor(plot_data$PairType, levels = c("P-P", "N-N", "N-P"))

# 比较组设定
comparisons <- list(c("P-P", "N-N"), c("P-P", "N-P"), c("N-N", "N-P"))

ggplot(plot_data, aes(x = PairType, y = Diff, fill = PairType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    tip.length = 0.01,      # 横杠高度（靠近星号）
    vjust = -0.3,           # 星号向上挪一点
    step.increase = 0.15    # 每一组比较整体往上抬高一点
  ) +
  scale_fill_manual(values = c("P-P" = "#66c2a5", "N-N" = "#8da0cb", "N-P" = "#fc8d62")) +
  labs(
    x = "Species Pair Type",
    y = "Breakpoints Time Difference"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

ggsave(filename = paste0("Breakpoints_Time_Violin_3Group.png"), width = 5, height = 4, dpi = 300)


# 
# ggplot(plot_data, aes(x = PairType, y = Diff, fill = PairType)) +
#   geom_violin(alpha = 0.6, trim = FALSE) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
#   scale_fill_manual(values = c("P-P" = "#66c2a5", "N-P" = "#fc8d62", "N-N" = "#8da0cb")) +
#   labs(x = "Pair Type", y = "Breakpoints Time Difference") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")
# 
# ggplot(plot_data, aes(x = PairType, y = Diff, fill = PairType)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
#   stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c("N-N", "P-P"), c("N-P", "P-P"), c("N-P", "N-N"))) +
#   scale_fill_manual(values = c("P-P" = "#66c2a5", "N-P" = "#fc8d62", "N-N" = "#8da0cb")) +
#   labs(x = "Pair Type", y = "Breakpoints Time Difference") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")
# 
# 
# 
# 
# library(ggplot2)
# library(ggpubr)
# 
# # 确保分组顺序
# plot_data$PairType <- factor(plot_data$PairType, levels = c("P-P", "N-N", "N-P"))
# 
# # 比较组设定
# comparisons <- list(c("P-P", "N-N"), c("P-P", "N-P"), c("N-N", "N-P"))
# 
# # 绘图
# ggplot(plot_data, aes(x = PairType, y = Diff, fill = PairType)) +
#   geom_violin(trim = FALSE, alpha = 0.7) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test",
#                      label = "p.signif", hide.ns = TRUE,
#                      tip.length = 0.03, vjust = -1) +  # 👈 向上移动星号
#   scale_fill_manual(values = c("P-P" = "#66c2a5", "N-N" = "#8da0cb", "N-P" = "#fc8d62")) +
#   labs(
#     x = "Species Pair Type",
#     y = "Breakpoints Time Difference"
#   ) +
#   theme_minimal(base_size = 15) +
#   theme(legend.position = "none")



