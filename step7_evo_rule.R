# Songyao 2024
#--------------------------------------原图------------------------------------------
# getwd() 
# rm(list = ls())
# library(ggsci) 
# library(ggplot2) 
# library(vioplot) 
# library(gridExtra)
# library(reshape2)
# # Set working directory
# dir="Y:/sda/results/Evolution_cortical_shape/statistic_landmarks/"
# setwd(dir)
# 
# # 定义向量
# primate <- c(31,11,7,8)
# non_pri <- c(10,1,11,10)

# # 保存
# png("Y:/sda/results/Evolution_cortical_shape/figure/pit_evo_rule.png", width = 700, height = 400, res = 100)
# # 设置绘图区域为1行2列
# par(mfrow = c(1, 2))  
# # 设置每个条的颜色
# color_palette <- pal_npg('nrc')(10)
# bar_colors <- color_palette[5:8]
# # 使用adjustcolor()函数设置透明度，alpha值为0.6 (60%透明度)
# bar_colors_alpha <- adjustcolor(bar_colors, alpha.f = 0.5)
# # 绘制第一个条形图 (Primate)
# barplot(primate, 
#         col = bar_colors_alpha,  # 设置颜色
#         names.arg = c("A≤Min","Max<A<Min", "A≥Max", "A=B1=B2"),  # X轴标签
#         main = "Primate",  # 图标题
#         ylab = "Count",  # Y轴标签
#         ylim = c(0, 35))  # 设置Y轴范围
# 
# # 绘制第二个条形图 (Non-Primate)
# barplot(non_pri, 
#         col = bar_colors_alpha,  # 设置颜色
#         names.arg = c("A≤Min","Max<A<Min", "A≥Max", "A=B1=B2"),  # X轴标签
#         main = "Non-Primate",  # 图标题
#         # ylab = "Count",  # Y轴标签
#         ylim = c(0, 35))  # 设置Y轴范围
# 
# dev.off()

#--------------------------------------2025-03-21------------------------------------------
# 清理环境
rm(list = ls())
library(ggplot2)
library(ggsci)
library(tidyr)   # 使用 pivot_longer()
library(dplyr)

# 设置工作目录
dir <- "Y:/sda/results/Evolution_cortical_shape/statistic_landmarks/"
setwd(dir)

# 定义数据
primate <- c(31, 11, 8, 7)
non_pri <- c(10, 1, 10, 11)

# 创建数据框
data <- data.frame(
  Category = c("A≤Min", "Max<A<Min", "A=B1=B2", "A≥Max"),
  Primate = primate,
  NonPrimate = non_pri
)

# 将数据转换为长格式
data_long <- pivot_longer(data, cols = -Category, 
                          names_to = "Group", values_to = "Count")

# 确保横轴顺序与数据一致
data_long$Category <- factor(data_long$Category, levels = data$Category)

# 配色方案
colors <- pal_npg("nrc")(10)  # Nature配色方案
primate_color <- adjustcolor(colors[4], alpha.f = 0.9)   # 深色 (Primate)
non_primate_color <- adjustcolor(colors[4], alpha.f = 0.4)  # 浅色 (Non-Primate)

# 绘图
p <- ggplot(data_long, aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Primate" = primate_color, "NonPrimate" = non_primate_color)) +
  labs(
    # title = "Distribution of Evolutionary Rules in Primate and Non-Primate",
    x = NULL,   # 不显示默认横轴标签
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),   # 横轴标签加粗
    axis.text.y = element_text(size = 14, face = "bold"),   # 纵轴刻度加粗
    axis.title.y = element_text(size = 16, face = "bold"),  # 纵轴标题加粗
    legend.position = c(0.95, 0.85),  # 图例在右上角
    legend.justification = c(1, 1),
    legend.title = element_blank()
  ) +
  # 在横轴中间添加类别标签
  annotate("text", x = 2.5, y = max(data_long$Count) + 2, 
           label = "Evolutionary Rule", fontface = "bold", size = 7)

# 保存图片
ggsave("Y:/sda/results/Evolution_cortical_shape/figure/pit_evo_rule_combined.png", 
       plot = p, width = 7, height = 4, dpi = 300)

# 显示图形
print(p)



# #--------------------------------------修改old------------------------------------------
# # 定义向量
# primate <- c(16, 38, 3) / 57
# non_pri <- c(2, 26, 4) / 32
# color_palette <- pal_npg('nrc')(10)
# 
# # 创建数据框
# data <- data.frame(
#   Category = c("<Min", "Midrange", ">Max" , "<Min", "Midrange", ">Max"),
#   Group = c("Primate", "Primate", "Primate", "Non-primate", "Non-primate", "Non-primate"),
#   Value = c(primate, non_pri)
# )
# 
# # 绘制分组柱状图
# p <- ggplot(data, aes(x = Category, y = Value, fill = Group)) +
#   geom_bar(stat = "identity", position = "dodge", alpha = 0.6) +
#   labs(x = "Category", y = "Value", fill = "Group") +
#   scale_fill_manual(values = c("Primate" = color_palette[1], "Non-primate" = color_palette[2])) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = "top"
#   )
# 
# 
# # 保存图形
# ggsave('Y:/sda/results/Evolution_cortical_shape/statistic_landmarks/pit_evo_rule_v2.png',
#        plot = p, width = 3, height = 3, dpi = 300)





























