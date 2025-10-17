# 清空环境
rm(list = ls())
setwd("Y:/results/Evolution_cortical_shape/statistic/kinship/")

library(R.matlab)
library(vegan)
library(ggplot2)
library(viridis)

# 读取数据
DIST_mat <- readMat("DIST_normalized.mat")$DIST.normalized  # 90x90 亲缘关系矩阵，已归一化到 [-1, 1]
# SIM_mat <- readMat("pit_corr_order.mat")$pit.corr.order     # 90x90 pit 相似性矩阵
SIM_mat <- readMat("pit_corr_Nring_order.mat")$pit.corr.order     # 90x90 pit 相似性矩阵
# 将相似性矩阵转为距离矩阵：距离 = 1 - 相似性
dist_phylo <- as.dist(1 - DIST_mat)
dist_pit <- as.dist(1 - SIM_mat)

# 计算 Mantel 检验
mantel_result <- mantel(dist_pit, dist_phylo, method = "pearson", permutations = 999)
print(mantel_result)

# 转换为散点图数据
get_upper_tri <- function(mat) {
  mat[lower.tri(mat, diag = TRUE)] <- NA
  return(mat)
}
DIST_mat_upper <- get_upper_tri(DIST_mat)
SIM_mat_upper <- get_upper_tri(SIM_mat)

df <- data.frame(
  kinship = as.vector(DIST_mat_upper),
  pit_sim = as.vector(SIM_mat_upper)
)
df <- na.omit(df)
# Mantel test 结果
mantel_r <- round(mantel_result$statistic, 3)
mantel_p <- mantel_result$signif

# 绘图
p <- ggplot(df, aes(x = kinship, y = pit_sim)) +
  geom_point(color = "#3C5488", size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "#E64B35", se = TRUE, linewidth = 1) +
  labs(
    x = "Phylogenetic Similarity",
    y = "Pit Distribution Similarity"
  ) +
  annotate("text",
           x = -0.9,
           y = 0.9,  # 固定在 y = 0.1
           hjust = 0,
           size = 5,
           color = "black",
           label = paste0("Mantel r = ", signif(mantel_r, 2), ", p = ", signif(mantel_p, 2))
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# 输出图形
print(p)

# 保存为 PNG 图片
ggsave(filename = "DIST_PitSim_mantel_plot.png", plot = p,
       width = 5, height = 4, dpi = 300)


# --------------------------------------------- 子模块 ---------------------------------------------
# 定义子模块编号
submodules <- list(1:45, 46:58, 61:71, 72:78, 79:86)
n_submodules <- length(submodules)

# 初始化 Mantel 检验结果矩阵
mantel_matrix <- matrix(NA, nrow = n_submodules, ncol = n_submodules)
rownames(mantel_matrix) <- colnames(mantel_matrix) <- paste0("M", 1:n_submodules)

# 计算子模块之间的 Mantel 检验
for (i in 1:(n_submodules - 1)) {
  for (j in (i + 1):n_submodules) {
    sub_dist_phylo <- as.dist(1 - DIST_mat[submodules[[i]], submodules[[j]]])
    sub_dist_pit <- as.dist(1 - SIM_mat[submodules[[i]], submodules[[j]]])
    sub_mantel_result <- mantel(sub_dist_pit, sub_dist_phylo, method = "pearson", permutations = 999)
    mantel_matrix[i, j] <- sub_mantel_result$statistic
  }
}

# 绘制热图（仅上三角矩阵）
get_upper_tri_heatmap <- function(mat) {
  mat[lower.tri(mat, diag = TRUE)] <- NA
  return(mat)
}
mantel_matrix_upper <- get_upper_tri_heatmap(mantel_matrix)

pheatmap(mantel_matrix_upper,
         color = viridis(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Mantel Test Results between Sub - modules",
         breaks = seq(min(mantel_matrix, na.rm = TRUE), max(mantel_matrix, na.rm = TRUE), length.out = 101))

# 保存热图为 PNG 图片
ggsave(filename = "submodules_mantel_heatmap.png",
       width = 5, height = 4, dpi = 300)





