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
library(grid)
library(R.matlab)
library(ggpubr)

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


#------------------------------内部相似性-------------------------------
#----------------------------24*24 vs 24*24-----------------------------
diff <- as.matrix(abun_unique)
same <- as.matrix(abun_shared)
cor_diff_t <- cor(diff, use = "pairwise.complete.obs", method = "pearson")
cor_same_t <- cor(same, use = "pairwise.complete.obs", method = "pearson")

heatmap_diff_t <- pheatmap(cor_diff_t, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("#7b3294", "#c2a5cf", "#f7f7f7", "#a6dba0", "#008837"))(100), 
                           show_rownames = TRUE,  # 显示基因名称
                           show_colnames = TRUE,  # 显示基因名称
                           labels_row = cell_types,  # 行标签使用基因名称
                           labels_col = cell_types,  # 列标签使用基因名称
                           main = "Specific Region")
heatmap_same_t <- pheatmap(cor_same_t, 
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE, 
                           color = colorRampPalette(c("#7b3294", "#c2a5cf", "#f7f7f7", "#a6dba0", "#008837"))(100), 
                           show_rownames = TRUE,
                           show_colnames = TRUE,  
                           labels_row = cell_types,  # 行标签使用基因名称
                           labels_col = cell_types,  # 列标签使用基因名称
                           main = "Shared Region")


combined_plot <- grid.arrange(heatmap_diff_t$gtable, heatmap_same_t$gtable, ncol = 2)
ggsave("combined_heatmap.png", plot = combined_plot, width = 10, height = 4.5, dpi = 300)


