# Songyao
# 小提琴图
# 加载必要的库
library(ggsci)
library(ggplot2)
library("R.matlab")
library(openxlsx)

# 设置输出路径
outdir = "Y:/sda/results/Evolution_cortical_shape/figure/"

# 读取数据
data <- read.xlsx('Y:/sda/results/Evolution_cortical_shape/statistic/saved/breakpoints_and_differentaiation_time.xlsx')
breakpoints <- data[[2]]

# 读取 DIST 矩阵
DIST <- readMat('Y:/sda/results/Evolution_cortical_shape/statistic/kinship/Species_kinship_DIST.mat')
DIST <- DIST[["DIST"]]

# 初始化变量
num_species <- length(breakpoints)
breakpoints_diff <- numeric()
dist_diff <- numeric()

# 计算两两物种间的差值和距离
for (i in 1:(num_species - 1)) {
  for (j in (i + 1):num_species) {
    breakpoints_diff <- c(breakpoints_diff, abs(breakpoints[i] - breakpoints[j]))
    dist_diff <- c(dist_diff, DIST[i, j])
  }
}

# 构建小提琴图数据框
violin_data <- data.frame(
  DIST = dist_diff,
  Breakpoint_Diff = breakpoints_diff
)
violin_data$DIST_Range <- cut(
  violin_data$DIST,
  breaks = c(0, 50, 100, 156),
  labels = c("0-50", "51-100", "101-156"),
  include.lowest = TRUE
)

# 计算中位数
median_values <- aggregate(Breakpoint_Diff ~ DIST_Range, data = violin_data, FUN = median)

# 绘制小提琴图
violin_plot <- ggplot(violin_data, aes(x = DIST_Range, y = Breakpoint_Diff, fill = DIST_Range, color = DIST_Range)) +
  geom_violin(alpha = 0.3, width = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black", alpha = 0.5) +
  geom_point(data = median_values, aes(x = DIST_Range, y = Breakpoint_Diff), 
             color = "white", size = 3, shape = 21, fill = "white", stroke = 1) +
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    x = "Distance Range",
    y = "Breakpoint Difference",
    fill = "DIST Range",
    color = "DIST Range"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, size = 0.4)
  )

# 打印图形
print(violin_plot)

# 保存小提琴图为PNG
ggsave(
  filename = paste0("DIST_vs_BreakpointsDiff_Violin.png"),
  plot = violin_plot,
  width = 5.5, height = 5.5, dpi = 600, bg = "white"
)


# # 加载必要的库
# library(ggsci)
# library(ggplot2)
# library("R.matlab")
# library(openxlsx)
# library(patchwork)  # 用于拼接图形
# 
# # 设置输出路径
# outdir = "Y:/sda/results/Evolution_cortical_shape/figure/"
# 
# # 读取数据
# data <- read.xlsx('Y:/sda/results/Evolution_cortical_shape/statistic/saved/breakpoints_and_differentaiation_time.xlsx')
# breakpoints <- data[[2]]
# 
# # 读取 DIST 矩阵
# DIST <- readMat('Y:/sda/results/Evolution_cortical_shape/statistic/kinship/Species_kinship_DIST.mat')
# DIST <- DIST[["DIST"]]
# 
# # 初始化变量
# num_species <- length(breakpoints)
# breakpoints_diff <- numeric()
# dist_diff <- numeric()
# 
# # 计算两两物种间的差值和距离
# for (i in 1:(num_species - 1)) {
#   for (j in (i + 1):num_species) {
#     breakpoints_diff <- c(breakpoints_diff, abs(breakpoints[i] - breakpoints[j]))
#     dist_diff <- c(dist_diff, DIST[i, j])
#   }
# }
# 
# # 整合为数据框
# plot_data_breakpoints <- data.frame(
#   Breakpoint_Diff = breakpoints_diff,
#   DIST = dist_diff
# )
# 
# # 将 DIST 划分为三组
# violin_data <- data.frame(
#   DIST = dist_diff,
#   Breakpoint_Diff = breakpoints_diff
# )
# violin_data$DIST_Range <- cut(
#   violin_data$DIST,
#   breaks = c(0, 50, 100, 156),
#   labels = c("0-50", "51-100", "101-156"),
#   include.lowest = TRUE
# )
# # 计算中位数并将其添加到数据框中
# median_values <- aggregate(Breakpoint_Diff ~ DIST_Range, data = violin_data, FUN = median)
# 
# # 绘制小提琴图
# violin_plot <- ggplot(violin_data, aes(x = DIST_Range, y = Breakpoint_Diff, fill = DIST_Range, color = DIST_Range)) +
#   geom_violin(alpha = 0.3, width = 0.7) +
#   geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black", alpha = 0.5) +
#   geom_point(data = median_values, aes(x = DIST_Range, y = Breakpoint_Diff), 
#              color = "white", size = 3, shape = 21, fill = "white", stroke = 1) +  # 用白色点表示中位数
#   scale_fill_npg() +
#   scale_fill_npg() +
#   scale_color_npg() +
#   labs(
#     x = "Distance Range",
#     y = "Breakpoint Difference",
#     fill = "DIST Range",
#     color = "DIST Range"
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(face = "bold", size = 12),
#     axis.title = element_text(face = "bold", size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "none",
#     plot.margin = margin(10, 10, 10, 10),
#     panel.border = element_rect(color = "black", fill = NA, size = 0.4),
#     geom_hline(yintercept = 0, linetype = "dashed", color = "black")
#   )
# 
# # 计算相关性r和线性拟合的统计量
# model <- lm(Breakpoint_Diff ~ DIST, data = plot_data_breakpoints)
# r_squared <- summary(model)$r.squared
# correlation <- cor(plot_data_breakpoints$DIST, plot_data_breakpoints$Breakpoint_Diff)
# p_value <- summary(model)$coefficients[2, 4]  # 提取斜率的p值
# 
# # 将R²、相关性r和p值转换为文本标签
# annotation_text <- paste(
#   "R² =", round(r_squared, 2), "\n",
#   "r =", round(correlation, 2), "\n",
#   "p =", format(p_value, scientific = TRUE, digits = 2)
# )
# 
# # 绘制线性拟合图并添加文本标签
# fit_plot <- ggplot(plot_data_breakpoints, aes(x = DIST, y = Breakpoint_Diff)) +
#   geom_point(color = "darkblue", alpha = 0.6, size = 0.3) +
#   geom_smooth(method = "lm", color = "red", se = TRUE) +
#   labs(
#     x = "Distance (DIST)",
#     y = "Breakpoint Time Difference"
#   ) +
#   annotate("text", x = max(plot_data_breakpoints$DIST) * 0.1, y = max(plot_data_breakpoints$Breakpoint_Diff) * 0.9,
#            label = annotation_text, color = "black", size = 6, hjust = 0) +
#   theme_minimal() +
#   theme(
#     axis.text = element_text(face = "bold", size = 12),
#     axis.title = element_text(face = "bold", size = 14),
#     plot.margin = margin(10, 10, 10, 10),
#     panel.border = element_rect(color = "black", fill = NA, size = 0.4)
#   )
# 
# # 将两个图拼接为一个图
# combined_plot <- fit_plot / violin_plot
# 
# # 保存拼接后的图形
# # ggsave(filename = paste0(outdir, "Violin_Fit_BREAKPOINT_vs_DIST.png"), 
#        #plot = combined_plot, width = 6, height = 10, dpi = 300)
# 
# # 打印拼接后的图形以供预览
# print(combined_plot)
