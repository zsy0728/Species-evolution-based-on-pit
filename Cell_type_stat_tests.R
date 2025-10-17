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
label <- as.vector(data_list$label.schaefer100) 
abun_full <- read.csv(file.path("Y:/results/Evolution_cortical_shape/data_info/schaeffer_Jorstad_100_7Net_expr_mat_new_NormZscore0.3.csv"))
abun <- t(as.matrix(abun_full[,1:50]))  # 50*24
unique_idx <- which(label < 0)
shared_idx <- which(label > 0)
# shared_idx <- setdiff(1:50, unique_idx)
abun_unique <- abun[unique_idx, ]   # 17*24
abun_shared <- abun[shared_idx, ]   # 19*24/33*24

# brain_regions <- abun_full[, 1:50]  
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
# ------------------------------------ 差值展示 ------------------------------------
unique_mean <- colMeans(abun_unique)
shared_mean <- colMeans(abun_shared)
abundance_diff <-unique_mean - shared_mean 
# percent_diff <- 100 * (unique_mean - shared_mean) / shared_mean
diff_data <- data.frame(
  cell_type = cell_types,
  percent_diff = abundance_diff,
  cell_group = cell_group
)
diff_data <- diff_data[order(-diff_data$percent_diff), ]
diff_data$cell_type <- factor(diff_data$cell_type, levels = diff_data$cell_type)

group_colors <- c(
  "interneurons" = "#FF77B7",
  "non-neuronal cells" = "#FFAF61",
  "excitatory neurons" = "#B1D690"
)
p <- ggplot(diff_data, aes(x = cell_type, y = percent_diff, fill = cell_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = group_colors, name = "Cell Group") +
  labs(
    title = "Difference in Cell Abundance",
    x = NULL
    #y = "Relative Difference (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# 7. 保存图像
ggsave("cell_abundance_diff_barplot.png", plot = p, width = 9, height = 4, dpi = 300)

# 8. 保存数据
#writexl::write_xlsx(diff_data, "cell_abundance_diff.xlsx")

# ------------------------------------ 统计检验 ------------------------------------
# ------------------------------------ 24 cell ------------------------------------
results_list <- list()

for (i in 1:ncol(abun)) {
  cell_name <- cell_types[i]
  group1 <- abun_unique[, i]
  group2 <- abun_shared[, i]
  
  # 计算均值和检验
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  test <- wilcox.test(group1, group2)
  
  results_list[[i]] <- data.frame(
    cell_type = cell_name,
    unique_mean = mean1,
    shared_mean = mean2,
    p_value = test$p.value
  )
}

# 合并为一个数据框
results_df <- do.call(rbind, results_list)
results_df$p_adj <- p.adjust(results_df$p_value, method = "fdr")

# 添加显著性标注
results_df$significance <- cut(results_df$p_adj,
                               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                               labels = c("***", "**", "*", "")
)

# 转长格式用于 ggplot 柱状图
library(tidyr)
plot_df <- results_df |>
  pivot_longer(cols = c("unique_mean", "shared_mean"),
               names_to = "region_type", values_to = "mean_abundance")
plot_df$region_type <- factor(plot_df$region_type, levels = c("shared_mean", "unique_mean"))

# 加入细胞分类信息
plot_df$cell_group <- rep(cell_group, each = 2)
plot_df$significance <- rep(results_df$significance, each = 2)

# 画图：按每个细胞类型作 unique / shared 的柱状图
library(ggplot2)
group_colors <- c(
  "interneurons" = "#FF77B7",
  "non-neuronal cells" = "#FFAF61",
  "excitatory neurons" = "#B1D690"
)

p <- ggplot(plot_df, aes(x = cell_type, y = mean_abundance, fill = region_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(
    data = subset(plot_df, region_type == "unique_mean"),
    aes(label = significance),
    vjust = -0.5, size = 5
  ) +
  facet_wrap(~cell_group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("shared_mean" = "#999999", "unique_mean" = "#377EB8"),
                    labels = c("Shared Region", "Unique Region"),
                    name = "Region Type") +
  labs(
    title = "Cell Type Abundance in Shared vs Unique Cortical Regions",
    x = NULL, y = "Mean Abundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title.y = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# 保存图像和数据
ggsave("cell_type_abundance_comparison.pdf", plot = p, width = 11, height = 4.5, dpi = 300)
writexl::write_xlsx(results_df, "cell_type_abundance_wilcox_results.xlsx")

# ------------------------------------ 3 Group ------------------------------------
group_labels <- unique(cell_group)
results <- data.frame()

for (grp in group_labels) {
  grp_idx <- which(cell_group == grp)
  
  unique_vals <- rowMeans(abun_unique[, grp_idx])
  shared_vals <- rowMeans(abun_shared[, grp_idx])
  
  # Wilcoxon检验
  p_val <- wilcox.test(unique_vals, shared_vals)$p.value
  sig <- ifelse(p_val < 0.001, "***",
                ifelse(p_val < 0.01, "**",
                       ifelse(p_val < 0.05, "*", "")))
  
  results <- rbind(results, data.frame(
    cell_group = grp,
    region_type = "Unique",
    mean_abundance = mean(unique_vals),
    p_value = p_val,
    significance = sig
  ))
  
  results <- rbind(results, data.frame(
    cell_group = grp,
    region_type = "Shared",
    mean_abundance = mean(shared_vals),
    p_value = p_val,
    significance = sig
  ))
}

# === 可视化 ===
results$cell_group <- factor(results$cell_group, levels = c("interneurons", "excitatory neurons", "non-neuronal cells"))
results$region_type <- factor(results$region_type, levels = c("Shared", "Unique"))

# 星号位置
star_df <- results %>%
  group_by(cell_group) %>%
  summarise(y_pos = max(mean_abundance) + 0.02,
            significance = first(significance))

p <- ggplot(results, aes(x = cell_group, y = mean_abundance, fill = region_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(data = star_df, aes(x = cell_group, y = y_pos, label = significance),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c("Shared" = "#999999", "Unique" = "#377EB8")) +
  labs(
    title = "Cell Class Abundance in Shared vs Unique Cortical Regions",
    x = NULL, y = "Mean Abundance", fill = "Region Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave("cell_group_abundance_comparison.pdf", plot = p, width = 8, height = 4, dpi = 300)

# ------------------------------------ 3 Group Percentage ------------------------------------
# ==== 将24种cell归为3类后求丰度 ====
group_labels <- c("interneurons", "excitatory neurons", "non-neuronal cells")
grouped_matrix <- matrix(0, nrow = nrow(abun), ncol = length(group_labels))
colnames(grouped_matrix) <- group_labels

for (g in group_labels) {
  grouped_matrix[, g] <- rowSums(abun[, cell_group == g])
}

# ==== 每行（每个脑区）做归一化：得到比例 ====
grouped_prop <- grouped_matrix / rowSums(grouped_matrix)

# ==== 加入区域标签 ====
region_type <- ifelse(label < 0, "Unique", ifelse(label > 0, "Shared", NA))
df_prop <- as.data.frame(grouped_prop)
df_prop$region_type <- region_type

# ==== 转换为长格式便于绘图 ====
df_long <- melt(df_prop, id.vars = "region_type", variable.name = "cell_class", value.name = "proportion")

# ==== 绘图 + 显著性统计 ====
# 统计检验
p_values <- df_long %>%
  group_by(cell_class) %>%
  summarise(
    p_value = wilcox.test(proportion ~ region_type)$p.value
  ) %>%
  mutate(
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 合并用于添加星号
df_long <- left_join(df_long, p_values, by = "cell_class")

# 确定每组上方标星位置
star_pos <- df_long %>%
  group_by(cell_class) %>%
  summarise(y_pos = max(proportion) + 0.03, signif = first(signif))

# ==== 画箱线图 ====
p <- ggplot(df_long, aes(x = cell_class, y = proportion, fill = region_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.6, size = 1) +
  geom_text(data = star_pos, aes(x = cell_class, y = y_pos, label = signif),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c("Shared" = "#4DAF4A", "Unique" = "#984EA3")) +
  labs(
    title = "Cell Class Proportion Across Cortical Region Types",
    x = NULL, y = "Relative Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

ggsave("cell_class_proportion_comparison.pdf", p, width = 6, height = 5, dpi = 300)



#---------------------一个柱状图-------------------------
# 计算均值
shared_means <- colMeans(abun_shared_agg)
unique_means <- colMeans(abun_unique_agg)

# 合并成数据框
df_prop <- data.frame(
  cell_group = names(shared_means),
  shared = shared_means,
  unique = unique_means
)

# 转为长格式 + 归一化为比例（加总为1）
library(tidyr)
df_long <- pivot_longer(df_prop, cols = c("shared", "unique"), names_to = "region_type", values_to = "abundance")

# 按region_type归一化比例
df_long <- df_long %>%
  group_by(region_type) %>%
  mutate(proportion = abundance / sum(abundance)) %>%
  ungroup()
# 设置颜色
group_colors <- c(
  "interneurons" = "#FF77B7",
  "non-neuronal cells" = "#FFAF61",
  "excitatory neurons" = "#B1D690"
)

# 画图
library(ggplot2)
p <- ggplot(df_long, aes(x = region_type, y = proportion, fill = cell_group)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = group_colors, name = "Cell Group") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of Cell Groups in Shared vs Unique Regions",
    x = NULL, y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12)
  )

print(p)
out_dir <- "Y:/results/Evolution_cortical_shape/homo_with_related_species/cell_type/"
ggsave(filename = paste0(out_dir, "cell_group_proportion_barplot.png"),
       plot = p, width = 6, height = 5, dpi = 300)

