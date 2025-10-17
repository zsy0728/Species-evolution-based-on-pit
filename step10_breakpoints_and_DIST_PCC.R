# Songyao 2025.04.17
# 适合多变量的分析，这个实验中只有两个变量
rm(list = ls())

library(ggsci) 
library(ggplot2) 
library(gridExtra)
library("R.matlab")
library(openxlsx)
library(linkET)
library(vegan)
library(ggtext)
library(cols4all)
library(tidyverse)

# Set working directory
outdir = "Y:/sda/results/Evolution_cortical_shape/figure/"

# ----------------------------------- load breakpoints_and_DIST-----------------------------------
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
# ----------------------------------- Prepare data for correlation and mantel test -----------------------------------
# 构建两个矩阵：一个为差异度量，一个为谱系距离
diff_df <- data.frame(breakpoints_diff = breakpoints_diff,
                      dist_diff = dist_diff)

# 计算 Pearson 相关性
cor2 <- correlate(diff_df)
corr2 <- cor2 %>% as_md_tbl()
write.csv(corr2, file = paste0( "pearson_correlate.csv"), row.names = TRUE)

# ----------------------------------- Mantel test -----------------------------------
# 为 Mantel test 创建距离矩阵（转为矩阵形式）
breakpoint_matrix <- matrix(0, nrow = num_species, ncol = num_species)
dist_matrix <- matrix(0, nrow = num_species, ncol = num_species)

for (i in 1:(num_species - 1)) {
  for (j in (i + 1):num_species) {
    val1 <- abs(breakpoints[i] - breakpoints[j])
    val2 <- DIST[i, j]
    breakpoint_matrix[i, j] <- val1
    breakpoint_matrix[j, i] <- val1
    dist_matrix[i, j] <- val2
    dist_matrix[j, i] <- val2
  }
}

# 转换为距离对象
breakpoint_dist <- as.dist(breakpoint_matrix)
dist_dist <- as.dist(dist_matrix)

# 运行 Mantel test
mantel_res <- vegan::mantel(breakpoint_dist, dist_dist, method = "pearson", permutations = 9999)
mantel_df <- data.frame(
  var1 = "breakpoints_diff",
  var2 = "dist_diff",
  r = mantel_res$statistic,
  p = mantel_res$signif
)

write.csv(mantel_df, file = paste0(outdir, "mantel_result.csv"), row.names = FALSE)

# 对 mantel 结果进行分组便于绘图
mantel2 <- mantel_df %>%
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                 labels = c("<0.25", "0.25-0.5", ">=0.5")),
         p = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                 labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">= 0.05")))

# ----------------------------------- 可视化 -----------------------------------
# 绘制热图
p4 <- qcorrplot(cor2,
                grid_col = "#00468BFF",
                "white", "#42B540FF",
                grid_size = 0.2,
                type = "upper",
                diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(colours = c("#00468BFF", "white", "#42B540FF"),
                       limits = c(-1, 1))

# 添加显著性标签
p5 <- p4 +
  geom_mark(size = 4,
            only_mark = TRUE,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05,
            colour = 'white')

# 添加 mantel 连线（这里只有一个连线，但为了统一格式保留）
p6 <- p5 +
  geom_couple(data = mantel2,
              aes(x = var1, xend = var2,
                  y = var1, yend = var2,
                  colour = p, size = r),
              curvature = nice_curvature())

# 美化连线
p7 <- p6 +
  scale_size_manual(values = c(1, 2, 3)) +
  scale_colour_manual(values = c4a('brewer.set2', 4)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantel's p",
                               override.aes = list(size = 5),
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(
    text = element_text(size = 16, family = "serif"),
    plot.title = element_text(size = 16, colour = "black", hjust = 0.5),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.text.y = element_text(size = 16, color = "black", vjust = 0.5, hjust = 1, angle = 0),
    axis.text.x = element_text(size = 16, color = "black", vjust = 0.5, hjust = 0.5, angle = 0)
  )

# 保存图像或直接展示
print(p7)
ggsave(paste0("breakpoints_and_DIST_pcc_mantel.png"), plot = p7, width = 8, height = 6)

