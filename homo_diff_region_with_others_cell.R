#-------------------------Cell Type-----------------------------
rm(list = ls())
dir <- "Y:/results/Evolution_cortical_shape/homo_with_related_species/cell_type/"
setwd(dir)

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(writexl)
library(pheatmap)
library(gridExtra)
library(R.matlab)

# 读取数据
data_list <- readMat("Y:/results/Evolution_cortical_shape/homo_with_related_species/homo_diff_region_label_schaefer100.mat")
names(data_list)
label <- as.vector(data_list$label.schaefer100) 
abun_full <- read.csv(file.path("Y:/results/Evolution_cortical_shape/data_info/schaeffer_Jorstad_100_7Net_expr_mat_new_NormZscore0.3.csv"))
abun <- t(as.matrix(abun_full[,1:50]))  
unique_idx <- which(label < 0)
shared_idx <- which(label > 0)
abun_unique <- abun[unique_idx, ]   
abun_shared <- abun[shared_idx, ]   

brain_regions <- abun_full[, 1:50]  
cell_types <- abun_full[, 101]       # 提取第 101 列，代表细胞类型
cell_group <- c("interneurons","interneurons",  "interneurons",  "interneurons", "interneurons",
                "excitatory neurons", "excitatory neurons", "excitatory neurons", "excitatory neurons",
                "non-neuronal cells", "non-neuronal cells", "interneurons","non-neuronal cells",
                "non-neuronal cells","non-neuronal cells","excitatory neurons","excitatory neurons",
                "excitatory neurons","excitatory neurons","excitatory neurons","interneurons",
                "interneurons","interneurons","interneurons")

# interneurons 粉色
# non-neuronal cells 橙色
# excitatory neurons 绿色
annotation_data <- data.frame(Group = cell_group)
rownames(annotation_data) <- cell_types 

diff <- abun_unique
same <- abun_shared

# 计算diff和same区域的平均丰度
diff_mean <- colMeans(diff) 
same_mean <- colMeans(same)



#------------------------------柱状图0714------------------------------
abundance_diff <- diff_mean - same_mean
colnames(diff) <- rownames(annotation_data)
# cell_types 是长度为 24 的字符向量，对应列名
cell_types <- colnames(diff)

# cell_group 是长度为 24 的分组标签向量，比如：
cell_group <- c(
  rep("excitatory neurons", 10),
  rep("interneurons", 8),
  rep("non-neuronal cells", 6)
)
# 构建数据框
diff_data <- data.frame(
  cell_type = cell_types,
  abundance_diff = abundance_diff,
  cell_group = cell_group
)

# 排序
diff_data <- diff_data[order(-diff_data$abundance_diff), ]
diff_data$cell_type <- factor(diff_data$cell_type, levels = diff_data$cell_type)

# 自定义颜色
group_colors <- c(
  "interneurons" = "#FF77B7",
  "non-neuronal cells" = "#FFAF61",
  "excitatory neurons" = "#B1D690"
)

# 画柱状图
library(ggplot2)
p <- ggplot(diff_data, aes(x = cell_type, y = abundance_diff, fill = cell_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = group_colors, name = "Cell Group") +
  labs(
    title = "Difference in Cell Abundance",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 19, hjust = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# 展示图
p

# 保存图
ggsave("cell_abundance_difference.png", plot = p, width = 8, height = 3, dpi = 300)

# 可选：保存数据
# install.packages("writexl") # 如未安装
library(writexl)
write_xlsx(diff_data, "diff_data.xlsx")



#------------------------------柱状图-------------------------------
#----------------------------pic 1 差值----------------------------
abundance_diff <- diff_mean - same_mean
# 创建数据框存储细胞类型、差值和分组信息
diff_data <- data.frame(
  cell_type = cell_types,
  abundance_diff = abundance_diff,
  cell_group = cell_group
)
# 按差值从大到小排序
diff_data <- diff_data[order(-diff_data$abundance_diff), ]

# 创建颜色映射
group_colors <- c(
  "interneurons" = "#FF77B7",
  "non-neuronal cells" = "#FFAF61",
  "excitatory neurons" = "#B1D690"
)

# 按细胞组给cell_type加标签
diff_data$cell_type <- factor(diff_data$cell_type, levels = diff_data$cell_type)

# 创建柱状图
p<-ggplot(diff_data, aes(x = cell_type, y = abundance_diff, fill = cell_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = group_colors, name = "Cell Group") +  # 设置图例标题
  labs(
    title = "Difference in Cell Abundance",
    x = NULL,
    y = NULL, # 隐藏纵轴标题
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),  # 放大纵坐标字体
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 14),  # 放大图例标题字体
    legend.text = element_text(size = 12),  # 放大图例字体
    plot.title = element_text(size = 19, hjust = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")
p

ggsave("cell_abundance_difference.png", plot = p, width = 8, height = 3, dpi = 300)
write_xlsx(diff_data, "diff_data.xlsx")

#------------------------------内部相似性-------------------------------

#----------------------------24*24 vs 24*24-----------------------------
diff <- as.matrix(diff)
same <- as.matrix(same)
cor_diff_t <- cor(t(diff), use = "pairwise.complete.obs", method = "pearson")
cor_same_t <- cor(t(same), use = "pairwise.complete.obs", method = "pearson")

heatmap_diff_t <- pheatmap(cor_diff_t, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("blue", "white", "red"))(100), 
                           show_rownames = TRUE,  # 显示基因名称
                           show_colnames = TRUE,  # 显示基因名称
                           labels_row = cell_types,  # 行标签使用基因名称
                           labels_col = cell_types,  # 列标签使用基因名称
                           main = "Similarity Matrix: Diff Region",
                           filename = "heatmap_cell_type_diff.png",  # 保存为正方形图像
                           height = 6, width = 6.5, dpi = 300)

heatmap_same_t <- pheatmap(cor_same_t, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("blue", "white", "red"))(100), 
                           show_rownames = TRUE,
                           show_colnames = TRUE,  
                           labels_row = cell_types,  # 行标签使用基因名称
                           labels_col = cell_types,  # 列标签使用基因名称
                           main = "Similarity Matrix: Non-Diff Region",
                           filename = "heatmap_cell_type_same.png",  # 保存为正方形图像
                           height = 6, width = 6.5, dpi = 300)


# 使用 grid.arrange 显示两个热图
combined_plot <- grid.arrange(heatmap_diff_t$gtable, heatmap_same_t$gtable, ncol = 2)

# 保存为 PNG 文件
ggsave("combined_heatmap.png", plot = combined_plot, width = 12, height = 6, dpi = 300)


#----------------------------24*24 vs 76*76-----------------------------
diff <- as.matrix(diff)
same <- as.matrix(same)
cor_diff <- cor(diff, use = "pairwise.complete.obs", method = "pearson")
cor_same <- cor(same, use = "pairwise.complete.obs", method = "pearson")


heatmap_diff_t <- pheatmap(cor_diff, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("blue", "white", "red"))(100), 
                           show_rownames = TRUE,  # 显示基因名称
                           show_colnames = TRUE,  # 显示基因名称
                           labels_row = colnames(diff),  # 行标签使用基因名称
                           labels_col = colnames(diff),  # 列标签使用基因名称
                           main = "Cell Type Similarity",
                           filename = "heatmap_region_diff_by_cell.png",  # 保存为正方形图像
                           height = 12, width = 13, dpi = 300)

heatmap_same_t <- pheatmap(cor_same, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("blue", "white", "red"))(100), 
                           show_rownames = TRUE,  # 显示基因名称
                           show_colnames = TRUE,  # 昲显示基因名称
                           labels_row = colnames(same_clean),  # 行标签使用基因名称
                           labels_col = colnames(same_clean),  # 列标签使用基因名称
                           main = "Cell Type Similarity",
                           filename = "heatmap_region_same_by_cell.png",  # 保存为正方形图像
                           height = 12, width = 13, dpi = 300)
#----------------------------原始丰度热图-----------------------------
heatmap_diff_t <- pheatmap(diff, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("blue", "white", "red"))(100), 
                           show_rownames = TRUE,  # 显示细胞类型名称
                           show_colnames = TRUE,  # 显示细胞类型名称
                           labels_row = cell_types,  # 行标签使用细胞类型名称
                           labels_col = colnames(diff),  # 列标签使用细胞类型名称
                           main = "Cell Type Abundance: Diff Region",
                           filename = "heatmap_diff_cell_type.png",  # 保存热图
                           height = 12, width = 13, dpi = 300,
                           breaks = seq(0, 0.08, length.out = 101))  # 设置颜色范围

# 创建并保存same的热图
heatmap_same_t <- pheatmap(same_clean, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("blue", "white", "red"))(100), 
                           show_rownames = TRUE,  # 显示细胞类型名称
                           show_colnames = TRUE,  # 显示细胞类型名称
                           labels_row = cell_types,  # 行标签使用细胞类型名称
                           labels_col = colnames(same_clean),  # 列标签使用细胞类型名称
                           main = "Cell Type Abundance: Non-Diff Region",
                           filename = "heatmap_same_cell_type.png",  # 保存热图
                           height = 12, width = 13, dpi = 300,
                           breaks = seq(0, 0.08, length.out = 101))





#----------------------------pic 1 比值----------------------------
# 计算 (diff - same) / same 的值
abundance_ratio <- (diff_mean - same_mean) / same_mean * 100

# 创建新数据框存储比值及分组信息
ratio_data <- data.frame(
  cell_type = cell_types,
  abundance_ratio = abundance_ratio,
  cell_group = cell_group
)

# 按比值从大到小排序
ratio_data <- ratio_data[order(-ratio_data$abundance_ratio), ]
ratio_data$cell_type <- factor(ratio_data$cell_type, levels = ratio_data$cell_type)

# 创建颜色映射
group_colors <- c(
  "interneurons" = "#FF77B7",
  "non-neuronal cells" = "#FFAF61",
  "excitatory neurons" = "#B1D690"
)

# 创建柱状图
p_ratio <- ggplot(ratio_data, aes(x = abundance_ratio, y = cell_type, fill = cell_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = group_colors) + #, guide = "none"
  scale_x_continuous(labels = scales::percent_format(scale = 1)) + # 转为百分比显示
  labs(
    x = "Relative Difference in Cell Abundance",
    y = NULL, # 隐藏纵轴标题
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

# 保存图片
ggsave("cell_abundance_ratio.png", plot = p_ratio, width = 5, height = 5, dpi = 300)

# ------------------------Mean------------------------------------------
# 创建一个数据框存储每种细胞类型的显著性检验结果
result <- data.frame(cell_type = cell_types, p_value = rep(NA, 24))

# 对每个细胞类型进行 Mann-Whitney U 检验
for (i in 1:24) {
  
  # 提取当前细胞类型在感兴趣脑区和非感兴趣脑区的丰度数据
  diff_region_abundance <- diff[i, ]
  same_region_abundance <- same[i, ]
  
  # 去掉 same_region_abundance 中的 NA 数据
  same_region_abundance <- same_region_abundance[!is.na(same_region_abundance)]
  
  # 将 diff_region_abundance 转换为向量
  diff_region_abundance <- unlist(diff_region_abundance)
  
  # 进行 Mann-Whitney U 检验
  # result$p_value[i] <- wilcox.test(diff_region_abundance, same_region_abundance)$p.value
  # 使用 t 检验来判断感兴趣区域和非感兴趣区域的差异
  # result$p_value[i] <- t.test(diff_region_abundance, same_region_abundance)$p.value
  # result$p_value[i] <- kruskal.test(diff_region_abundance ~ same_region_abundance)$p.value
  result$p_value <- p.adjust(result$p_value, method = "fdr")  # 使用 FDR 校正
}

# 查看显著性检验结果
print(result)

# 对 p 值进行多重比较校正
result$adj_p_value <- p.adjust(result$p_value, method = "BH")

# 查看校正后的结果
print(result)

# -------------------------display mean-----------------------------------
mean_diff_regions <- rowMeans(brain_regions[, diff_region], na.rm = TRUE)
other_regions <- brain_regions[, -diff_region]
mean_other_regions <- rowMeans(other_regions, na.rm = TRUE)

# 创建一个数据框
mean_matrix <- data.frame(
  "Diff Regions Mean" = mean_diff_regions,
  "Other Regions Mean" = mean_other_regions
)

mean_matrix_melted <- melt(mean_matrix)
mean_matrix_melted$Cell_Type <- rep(abun[, 101], times = ncol(mean_matrix))  # 修正为times
mean_matrix_melted$Cell_Type <- factor(mean_matrix_melted$Cell_Type, levels = unique(abun[, 101]))

# 绘制热图
ggplot(mean_matrix_melted, aes(x = variable, y = Cell_Type, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#B2182B", midpoint = 0) +  # 使用 RdBu 配色方案
  theme_minimal() +
  xlab("Regions") +  # x轴标签
  ylab("Cell Types") +  # y轴标签
  theme(axis.text.x = element_text(angle = 90, hjust = 1), # 旋转x轴标签
        axis.text.y = element_text(size = 8)) +  # y轴标签大小
  scale_y_discrete(limits = rev(levels(mean_matrix_melted$Cell_Type))) +  # 确保y轴按照Cell_Type顺序
  theme(axis.text.y = element_text(size = 8))  # 控制纵轴字体大小

#------------------------------E:I在diff和non-diff之间是否有差异-------------------------------
# 获取细胞类型的分组索引
excitatory_neurons_idx <- which(cell_group == "excitatory neurons")
interneurons_idx <- which(cell_group == "interneurons")

# diff 数据集统计
diff_excitatory_sum <- rowSums(diff[excitatory_neurons_idx, ], na.rm = TRUE)
diff_interneurons_sum <- rowSums(diff[interneurons_idx, ], na.rm = TRUE)
diff_ratio <- sum(diff_excitatory_sum) / sum(diff_interneurons_sum)

# same_clean 数据集统计
same_clean_excitatory_sum <- rowSums(same_clean[excitatory_neurons_idx, ], na.rm = TRUE)
same_clean_interneurons_sum <- rowSums(same_clean[interneurons_idx, ], na.rm = TRUE)
same_clean_ratio <- sum(same_clean_excitatory_sum) / sum(same_clean_interneurons_sum)

# 输出比值
cat("diff 中 excitatory neurons / interneurons 的比值:", diff_ratio, "\n")
cat("same_clean 中 excitatory neurons / interneurons 的比值:", same_clean_ratio, "\n")
























