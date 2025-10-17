# 新版
rm(list = ls())
graphics.off()  # 关闭所有打开的图形设备
setwd("Y:/results/Evolution_cortical_shape/statistic/kinship/")
library(R.matlab)
library(readxl)
library(pheatmap)

# 加载数据
SIM_mat <- readMat("pit_corr_Nring_order.mat")$pit.corr.order     # 90x90 pit 相似性矩阵
# SIM_mat <- readMat("pit_corr_order.mat")$pit.corr.order
DIST_mat <- readMat("DIST_normalized.mat")$DIST.normalized
order_idx <- as.vector(readMat("phytree_indices.mat")$phytree.indices)

infodir = 'Y:/results/Evolution_cortical_shape/data_info/'
species_info <- readxl::read_excel(file.path(infodir, "Species_info.xlsx"))
species_info <- species_info[order_idx, ]
species_names_ordered <- species_info$Species

# 添加行列名（不调整矩阵顺序）
rownames(SIM_mat) <- colnames(SIM_mat) <- species_names_ordered
rownames(DIST_mat) <- colnames(DIST_mat) <- species_names_ordered

# 保存 PIT 相似性热图（默认配色）
# png("PIT_similarity_ring_heatmap_default.png", width = 2400, height = 2400, res = 300)
pdf("PIT_similarity_ring_heatmap_default.pdf", width = 8, height = 8)
pheatmap(SIM_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         main = "PIT Similarity Matrix")
dev.off()

# 保存 Phylogenetic 相似性热图（默认配色）
# png("Phylogenetic_ring_similarity_heatmap_default.png", width = 2400, height = 2400, res = 300)
pdf("Phylogenetic_ring_similarity_heatmap_default.pdf", width = 8, height = 8)
pheatmap(DIST_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         main = "Phylogenetic Similarity Matrix")
dev.off()



# ---------------------------------- 为了画图的heatmap --------------------------------
rm(list = ls())
graphics.off()  # 关闭所有打开的图形设备
setwd("Y:/sda/results/Evolution_cortical_shape/statistic/kinship/")
library(R.matlab)
library(readxl)
library(pheatmap)

# 加载数据
SIM_mat <- readMat("pit_corr_Nring_order.mat")$pit.corr.order     # 90x90 pit 相似性矩阵
# SIM_mat <- readMat("pit_corr_order.mat")$pit.corr.order
DIST_mat <- readMat("DIST_normalized.mat")$DIST.normalized
order_idx <- as.vector(readMat("phytree_indices.mat")$phytree.indices)

infodir = 'Y:/sda/results/Evolution_cortical_shape/data_info/'
species_info <- readxl::read_excel(file.path(infodir, "Species_info.xlsx"))
species_info <- species_info[order_idx, ]
species_names_ordered <- species_info$Species

# 添加行列名（不调整矩阵顺序）
rownames(SIM_mat) <- colnames(SIM_mat) <- species_names_ordered
rownames(DIST_mat) <- colnames(DIST_mat) <- species_names_ordered

# 定义权重系数，alpha 越小，new_SIM_mat 越接近 DIST_mat
alpha <- 0.9

# 叠加 DIST_mat 到 SIM_mat
new_SIM_mat <- alpha * SIM_mat + (1 - alpha) * DIST_mat

# 归一化 new_SIM_mat 到 0 到 1 之间
min_value <- min(new_SIM_mat)
max_value <- max(new_SIM_mat)
normalized_SIM_mat <- (new_SIM_mat - min_value) / (max_value - min_value)

# 保存归一化后的 SIM_mat 热图（默认配色）
pdf("Normalized_PIT_similarity_ring_heatmap.pdf", width = 8, height = 8)
pheatmap(normalized_SIM_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         main = "Normalized PIT Similarity Matrix")
dev.off()