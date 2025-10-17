# Songyao 2025.04.17 - DIST vs Breakpoints Difference 散点图 + 分三组拟合结果。用不上
rm(list = ls())

library(ggsci)
library(ggplot2)
library(R.matlab)
library(openxlsx)
library(vegan)
library(dplyr)

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

# ------------------ Compute correlation stats per group ------------------
annot_df <- plot_data %>%
  group_by(PairType) %>%
  summarise(
    r = cor(Distance, Diff, method = "pearson"),
    p = cor.test(Distance, Diff)$p.value,
    R2 = summary(lm(Diff ~ Distance))$r.squared
  ) %>%
  mutate(
    label = paste0("r = ", signif(r, 2),
                   "\np = ", signif(p, 2),
                   "\nR² = ", signif(R2, 2))
  )

# ------------------ Plot ------------------
p <- ggplot(plot_data, aes(x = Distance, y = Diff)) +
  geom_point(aes(color = PairType), size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  facet_wrap(~PairType, scales = "free") +
  geom_text(
    data = annot_df,
    mapping = aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.1, inherit.aes = FALSE,
    size = 4.5
  ) +
  scale_color_manual(values = c("P-P" = "#66c2a5", "N-P" = "#fc8d62", "N-N" = "#8da0cb")) +
  labs(
    x = "Evolutionary Distance (DIST)",
    y = "Breakpoints Time Difference",
    color = "Pair Type"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

print(p)

# 保存图像
# ggsave(filename = paste0(outdir, "DIST_vs_BreakpointsDiff_byGroup_annotated.png"),
#        plot = p, width = 10, height = 4, dpi = 300)
