# ------------------------------- 只画层次聚类结果 -------------------------------
rm(list = ls())
setwd("Y:/sda/results/Evolution_cortical_shape/statistic/kinship/")

library(R.matlab)
library(readxl)
library(pheatmap)

# 加载数据
SIM_mat <- readMat("pit_corr_Nring.mat")$distance.ring    
infodir <- 'Y:/sda/results/Evolution_cortical_shape/data_info/'
order_idx <- readMat("phytree_indices.mat")$phytree.indices |> as.vector()
species_info <- read_excel(file.path(infodir, "Species_info.xlsx"))
species_names <- species_info$Species

# 添加行列名（不调整矩阵顺序）
rownames(SIM_mat) <- colnames(SIM_mat) <- species_names

# 进行层次聚类
d <- dist(SIM_mat)
hc <- hclust(d, method = "complete")

# 绘制热图
pdf("Hierarchical Clustering for the pit Distribution.pdf", width = 12, height = 12)
pheatmap(SIM_mat, 
         cluster_rows = hc, 
         cluster_cols = hc,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Hierarchical Clustering for the pit Distribution",
         fontsize_row = 6,  # 可调整行名的字体大小
         fontsize_col = 6,   # 可调整列名的字体大小
         angle_col = 90
)
dev.off()

# ------------------------------- 层次聚类后色块的结果 -------------------------------
# 加载必要包
library(R.matlab)
library(readxl)
library(tidyverse)
library(cowplot)
library(pheatmap)

# 加载数据
infodir <- 'Y:/sda/results/Evolution_cortical_shape/data_info/'
name_list <- read.table(file.path(infodir, "phyTree_name_list.txt"), header = FALSE, stringsAsFactors = FALSE)
node_ages <- read_excel(file.path(infodir, "node_ages.xlsx"))
species_info <- read_excel(file.path(infodir, "Species_info.xlsx"))

# 提取变量
species_names <- species_info$Species
social_style <- as.character(species_info$Code_SocialStyle)
order_group <- ifelse(species_info$Order == "Primata", "Primata", "Other")
# 提取三个新的行为变量
fossorial <- as.character(species_info$Code_Fossorial)
terrestrial <- as.character(species_info$Code_Terrestrial)
arboreal <- as.character(species_info$Code_Arboreal)

# 主颜色设置
social_main <- "#f0a0e9"     # 深粉
fossorial_main <- "#2d8ba8"  # 深蓝
terrestrial_main <- "#6b9b72" # 深绿
arboreal_main <- "#b39cd0"   # 深橙
# 设置颜色函数（透明度处理）
apply_alpha <- function(hex, alpha) {
  rgb_val <- col2rgb(hex) / 255
  rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = alpha)
}
# 转成颜色向量
social_col <- ifelse(social_style == "1", social_main, apply_alpha(social_main, 0.5))
fossorial_col <- ifelse(fossorial == "1", fossorial_main, apply_alpha(fossorial_main, 0.5))
terrestrial_col <- ifelse(terrestrial == "1", terrestrial_main, apply_alpha(terrestrial_main, 0.5))
arboreal_col <- ifelse(arboreal == "1", arboreal_main, apply_alpha(arboreal_main, 0.5))
order_col <- c("Primata" = "#E64B35FF", "Other" = "#4DBBD5FF")  # 分类群：红 / 蓝

# 统一创建5个变量的长表格
plot_df <- rbind(
  tibble(Var = "Social Style",    y = 2, Index = 1:90, FillColor = social_col),
  tibble(Var = "Order",           y = 1, Index = 1:90, FillColor = order_col[order_group]),
  tibble(Var = "Fossorial",       y = 3, Index = 1:90, FillColor = fossorial_col),
  tibble(Var = "Terrestrial",     y = 4, Index = 1:90, FillColor = terrestrial_col),
  tibble(Var = "Arboreal",        y = 5, Index = 1:90, FillColor = arboreal_col)
)

# 进行层次聚类
# 这里假设存在一个距离矩阵 SIM_mat 用于聚类，如果没有需要根据实际情况构建
# 示例：假设已经有一个合适的 SIM_mat 矩阵
SIM_mat <- readMat("pit_corr_Nring.mat")$distance.ring    
rownames(SIM_mat) <- colnames(SIM_mat) <- species_names
d <- dist(SIM_mat)
hc <- hclust(d, method = "complete")

# 根据聚类结果对物种进行重新排序
sorted_indices <- hc$order
species_info <- species_info[sorted_indices, ]
species_names <- species_info$Species

# 根据新的排序更新 plot_df 中的 Index 和物种相关信息
plot_df$Index <- match(plot_df$Index, sorted_indices)
plot_df <- plot_df[order(plot_df$Index), ]

# 主图绘制
p_main <- ggplot(plot_df, aes(x = Index, y = y)) +
  geom_tile(aes(fill = FillColor), color = "white", width = 1, height = 1) +
  
  # 名字位置调整
  geom_text(
    data = tibble(x = 1:90, label = species_names),
    aes(x = x, y = 0.3, label = label),
    angle = 90, hjust = 1, size = 4,
    inherit.aes = FALSE
  ) +
  
  scale_fill_identity() +
  theme_void() +
  
  coord_cartesian(ylim = c(-11, 8)) +  # 修改范围包括5行色块 + 物种名空隙
  theme(plot.margin = margin(10, 30, 10, 30))

print(p_main)
ggsave("species_blocks_plot_sorted.pdf", p_main, width = 20, height = 4)

#  ------------------------------- 提取只包含 Order 的数据 ------------------------------- 

# 创建一个数据框用于绘图
df <- data.frame(
  x = seq_along(species_names),
  y = 0.3,  # 与原主图保持一致
  label = species_names
)

# 绘制旋转 90° 的物种名
p <- ggplot(df, aes(x = x, y = y, label = label)) +
  geom_text(
    angle = 90, 
    hjust = 0,  # 修改为左侧对齐
    size = 4,   # 与原主图保持一致
    vjust = 0.5 # 垂直居中
  ) +
  theme_void() +
  theme(plot.margin = margin(10, 30, 10, 30)) +  # 与原主图保持一致
  coord_cartesian(clip = "off")

# 保存为 PNG 文件
ggsave("hierarchical_clustered_species_names.png", p, width = 20, height = 4.5, dpi = 300)

#  ------------------------------- 只保存物种名顺序 ------------------------------- 
# 创建一个数据框用于绘图
df <- data.frame(
  x = seq_along(species_names),
  y = 1,
  label = species_names
)

# 绘制旋转 90° 的物种名
p <- ggplot(df, aes(x = x, y = y, label = label)) +
  geom_text(angle = 90, hjust = 0, vjust = 0.5, size = 3) +
  theme_void() +
  theme(plot.margin = margin(10, 30, 10, 30)) +
  coord_cartesian(clip = "off")

# 保存为 PNG 文件
ggsave("hierarchical_clustered_species_names.png", p, width = 15, height = 3, dpi = 300)



# ------------------------------- 衡量层次聚类结果 -------------------------------
# library(mclust)
# 
# # 自定义函数计算兰德指数
# rand_index <- function(cluster1, cluster2) {
#   n <- length(cluster1)
#   a <- 0
#   d <- 0
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       same_cluster1 <- cluster1[i] == cluster1[j]
#       same_cluster2 <- cluster2[i] == cluster2[j]
#       if (same_cluster1 & same_cluster2) {
#         a <- a + 1
#       } else if (!same_cluster1 & !same_cluster2) {
#         d <- d + 1
#       }
#     }
#   }
#   total_pairs <- choose(n, 2)
#   ri <- (a + d) / total_pairs
#   return(ri)
# }
# 
# # 示例聚类结果
# cluster1 <- c(0, 0, 1, 1, 2, 2)
# cluster2 <- c(0, 0, 1, 1, 1, 2)
# 
# # 计算兰德指数
# ri <- rand_index(cluster1, cluster2)
# cat("兰德指数 (RI):", ri, "\n")
# 
# # 计算调整兰德指数
# ari <- adjustedRandIndex(cluster1, cluster2)
# cat("调整兰德指数 (ARI):", ari, "\n")





# # 生成示例数据
# # 假设聚类结果，共有3个簇，10个样本
# cluster_result <- c(0, 1, 1, 1, 1, 0, 0, 0, 1, 1)
# # 假设真实标签，共有2个类别
# true_label <- c(0, 0, 1, 1, 1, 0, 0, 0, 1, 1)
# 
# # 计算聚类熵的函数
# calculate_cluster_entropy <- function(cluster_result, true_label) {
#   unique_clusters <- unique(cluster_result)
#   n <- length(true_label)
#   total_entropy <- 0
#   
#   for (cluster in unique_clusters) {
#     cluster_indices <- which(cluster_result == cluster)
#     n_i <- length(cluster_indices)
#     label_counts <- table(true_label[cluster_indices])
#     p_ij <- label_counts / n_i
#     
#     cluster_entropy <- -sum(p_ij * log2(p_ij))
#     total_entropy <- total_entropy + (n_i / n) * cluster_entropy
#   }
#   
#   return(total_entropy)
# }
# 
# # 调用函数计算聚类熵
# entropy <- calculate_cluster_entropy(cluster_result, true_label)
# cat("聚类熵为:", entropy, "\n")
# 





