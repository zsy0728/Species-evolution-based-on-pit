# 清理环境
rm(list = ls())
library(ggsci)
library(ggplot2)
library(gridExtra)
library("R.matlab")
library(openxlsx)
library(mgcv)     # GAM模型
library(segmented) # 分段回归模型

# 设置工作目录
outdir <- "Y:/sda/results/Evolution_cortical_shape/figure/"

# 加载数据
sp <- read.xlsx('Y:/sda/results/Evolution_cortical_shape/data_info/Species_info.xlsx')
na <- read.xlsx('Y:/sda/results/Evolution_cortical_shape/data_info/node_ages.xlsx')
pit_num <- readMat('Y:/sda/results/Evolution_cortical_shape/statistic_landmarks/number_of_pits_species_LH.mat')
all_paths <- readMat('Y:/sda/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
ages <- as.numeric(na[, 3])

# 将 sp 第 4 列中非 "Primata" 标记为 "Non-Primata"
sp[, 4] <- ifelse(sp[, 4] == "Primata", "Primata", "Non-Primata")

# 合并所有物种的数据
combined_primata <- data.frame(Age = numeric(), Count = numeric(), Species = character())
combined_non_primata <- data.frame(Age = numeric(), Count = numeric(), Species = character())

# 存储所有物种的分段拟合线数据
segmented_lines <- list()

# 合并数据 + 分段拟合
for (path_num in 1:90) {
  path <- all_paths[["all.paths"]][[path_num]][[1]]
  age <- ages[path]
  counts <- numeric()
  
  if (length(age) < 5) next  # 跳过数据点不足的情况
  
  for (j in 1:length(age)) {
    counts[j] <- pit_num[["pit.num.LH"]][path[j]]
  }
  
  # 数据归一化到 [0, 1] 区间
  counts_norm <- (counts - min(counts)) / (max(counts) - min(counts))
  # counts_norm <- (counts - 4) / (max(counts) - 4)
  
  species_name <- sp[path_num, 1]   # 物种名称
  group <- sp[path_num, 4]          # 分组
  is_primata <- group == "Primata"
  
  # 合并数据
  data <- data.frame(Age = age, Count = counts_norm, Species = species_name)
  
  # 分段回归拟合
  lm_model <- lm(Count ~ Age, data = data)
  seg_model <- segmented(lm_model, seg.Z = ~Age, npsi = 1, control = seg.control(quant = TRUE))
  
  # 存储分段拟合曲线数据
  seg_line <- data.frame(Age = data$Age, Fitted = fitted(seg_model), Group = group, Species = species_name)
  segmented_lines[[length(segmented_lines) + 1]] <- seg_line
  
  if (is_primata) {
    combined_primata <- rbind(combined_primata, data)
  } else {
    combined_non_primata <- rbind(combined_non_primata, data)
  }
}

# 将分段拟合线合并为数据框
segmented_lines_df <- do.call(rbind, segmented_lines)

# 对群组数据进行 GAM 拟合
gam_primata <- gam(Count ~ s(Age, k = 3), data = combined_primata)
gam_non_primata <- gam(Count ~ s(Age, k = 3), data = combined_non_primata)

# 预测拟合值
combined_primata$Fitted <- predict(gam_primata)
combined_non_primata$Fitted <- predict(gam_non_primata)

# 计算 R²
r2_primata <- summary(gam_primata)$r.sq
r2_non_primata <- summary(gam_non_primata)$r.sq

# 配色方案
palette <- pal_npg()(4)

# 预测拟合值和95%置信区间
pred_primata <- predict(gam_primata, se.fit = TRUE)
combined_primata$Fitted <- pred_primata$fit
combined_primata$Upper <- pred_primata$fit + 1.96 * pred_primata$se.fit
combined_primata$Lower <- pred_primata$fit - 1.96 * pred_primata$se.fit

pred_non_primata <- predict(gam_non_primata, se.fit = TRUE)
combined_non_primata$Fitted <- pred_non_primata$fit
combined_non_primata$Upper <- pred_non_primata$fit + 1.96 * pred_non_primata$se.fit
combined_non_primata$Lower <- pred_non_primata$fit - 1.96 * pred_non_primata$se.fit

# 绘制带95%置信区间的群组拟合图
p_primata <- ggplot(combined_primata, aes(x = Age, y = Count)) +
  geom_point(color = "black", size = 1, alpha = 0.3) +  # 散点图
  geom_line(data = segmented_lines_df[segmented_lines_df$Group == "Primata", ],
            aes(x = Age, y = Fitted, group = Species), 
            color = "grey70", size = 0.5, alpha = 0.6) +  # 分段回归线
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = palette[1], alpha = 0.4) +  # 95%置信区间
  geom_line(aes(y = Fitted), color = palette[1], size = 1.5) +  # GAM拟合曲线
  scale_x_reverse(limits = c(80, 0)) +  # X轴反转
  scale_y_continuous(limits = c(0, 1)) +  # Y轴归一化
  ylab("Normalized Count") +  # 修改纵坐标标题
  theme_minimal() +
  ggtitle("Primata") +
  theme(
    plot.title = element_text(color = "black", hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  annotate("text", x = 60, y = 0.9, 
           label = paste("R² =", round(r2_primata, 2)), 
           color = palette[4], size = 7, fontface = "bold")

p_non_primata <- ggplot(combined_non_primata, aes(x = Age, y = Count)) +
  geom_point(color = "black", size = 1, alpha = 0.3) +  # 散点图
  geom_line(data = segmented_lines_df[segmented_lines_df$Group == "Non-Primata", ],
            aes(x = Age, y = Fitted, group = Species), 
            color = "grey70", size = 0.5, alpha = 0.6) +  # 分段回归线
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = palette[2], alpha = 0.4) +  # 95%置信区间
  geom_line(aes(y = Fitted), color = palette[2], size = 1.5) +  # GAM拟合曲线
  scale_x_reverse(limits = c(80, 0)) +  # X轴反转
  scale_y_continuous(limits = c(0, 1)) +  # Y轴归一化
  ylab("Normalized Count") +  # 修改纵坐标标题
  theme_minimal() +
  ggtitle("Non-Primate") +
  theme(
    plot.title = element_text(color = "black", hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  annotate("text", x = 60, y = 0.9, 
           label = paste("R² =", round(r2_non_primata, 2)), 
           color = palette[4], size = 7, fontface = "bold")

# 将两组图上下排列
combined_plot <- grid.arrange(p_primata, p_non_primata, nrow = 2)

# 保存图片
ggsave(filename = paste0(outdir, "pit_num_GAM_group.png"), 
       plot = combined_plot, 
       width = 5, height = 8, dpi = 300)

# 将两组图排列在一起
combined_plot <- grid.arrange(p_primata, p_non_primata, nrow = 2)

# 保存图片
ggsave(filename = paste0(outdir, "pit_num_GAM_group_with_segmented_lines.png"), 
       plot = combined_plot, 
       width = 5, height = 8, dpi = 300)

