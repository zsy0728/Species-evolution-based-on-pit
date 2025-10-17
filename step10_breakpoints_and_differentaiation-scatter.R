# Songyao 2025.04.17 - DIST vs Breakpoints Difference 散点图
rm(list = ls())

library(ggsci)
library(ggplot2)
library(R.matlab)
library(openxlsx)
library(vegan)

# Set working directory
outdir = "Y:/sda/results/Evolution_cortical_shape/figure/"

# ------------------ Load data ------------------
data <- read.xlsx('Y:/sda/results/Evolution_cortical_shape/statistic/saved/breakpoints_and_differentaiation_time.xlsx')
breakpoints <- data[[2]]

DIST <- readMat('Y:/sda/results/Evolution_cortical_shape/statistic/kinship/Species_kinship_DIST.mat')
DIST <- DIST[["DIST"]]

num_species <- 90
breakpoints_diff <- numeric()
dist_diff <- numeric()

# 计算两两物种间的差值和距离
for (i in 1:(num_species - 1)) {
  for (j in (i + 1):num_species) {
    breakpoints_diff <- c(breakpoints_diff, abs(breakpoints[i] - breakpoints[j]))
    dist_diff <- c(dist_diff, DIST[i, j])
  }
}

# 整合为数据框
plot_data <- data.frame(
  Diff = breakpoints_diff,
  Distance = dist_diff
)

# ------------------ Plot ------------------
# 线性拟合模型
lm_model <- lm(Diff ~ Distance, data = plot_data)
lm_summary <- summary(lm_model)

# 拟合优度 R²
r_squared <- lm_summary$r.squared

# Pearson 相关性和 p 值
cor_result <- cor.test(plot_data$Diff, plot_data$Distance, method = "pearson")
r_value <- cor_result$estimate
p_value <- cor_result$p.value

# 绘图
p <- ggplot(plot_data, aes(x = Distance, y = Diff)) +
  geom_point(color = "#3C5488", size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "#E64B35", se = TRUE, linewidth = 1) +
  labs(
    x = "Evolutionary Distance (DIST)",
    y = "Breakpoints Time Difference"
  ) +
  annotate("text",
           x = max(plot_data$Distance) * 0.05,
           y = max(plot_data$Diff) * 0.85,
           hjust = 0,
           size = 5,
           color = "black",
           label = paste0("r = ", signif(r_value, 2), 
                          "\np = ", signif(p_value, 2),
                          "\nR² = ", signif(r_squared, 2))
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# 保存图像
#ggsave(filename = paste0("DIST_vs_BreakpointsDiff.png"),
   #    plot = p, width = 5, height = 4, dpi = 300)
