# Pit number increase vs. surface area expansion
# 
#------------------------------------------------------
rm(list = ls())
dir <- "Y:/sda/results/Evolution_cortical_shape/statistic/evo_speedy/"
setwd(dir)

library(ggplot2)
library(RColorBrewer)
library("R.matlab")

pri <- readMat(file.path("pitnum_growth_coupled_SA_expansion_primate_lh.mat"))
non <- readMat(file.path("pitnum_growth_coupled_SA_expansion_nonprimate_lh.mat"))
pri_field <- names(pri)[1]  # 提取灵长类字段名
non_field <- names(non)[1]  # 提取非灵长类字段名
pri_vec <- as.vector(pri[[pri_field]])
non_vec <- as.vector(non[[non_field]])

region_names <- c("Visual", "SMN", "DAN", "VAN", "Limbic", "FPN", "DMN")

# 创建数据框
data <- data.frame(
  Region = factor(rep(region_names, each = 2), levels = rev(region_names)),  # 保持脑区顺序
  Value = c(rbind(pri_vec, non_vec)),   
  Group = rep(c("Primate", "NonPrimate"), times = length(region_names))
)

# RGB 转 HEX 转换函数
rgb_to_hex <- function(r, g, b) {
  rgb(r / 255, g / 255, b / 255)
}
# 将 RGB 转换为 HEX
region_colors <- c(
  rgb_to_hex(120, 18, 134),   # Visual
  rgb_to_hex(70, 130, 180),   # SMN
  rgb_to_hex(0, 118, 14),     # DAN
  rgb_to_hex(196, 58, 250),   # VAN
  rgb_to_hex(220, 248, 164),  # Limbic
  rgb_to_hex(230, 148, 34),   # FPN
  rgb_to_hex(205, 62, 78)     # DMN
)

# 生成灵长类和非灵长类颜色
colors_primate <- region_colors          # 灵长类使用原色
colors_nonprimate <- adjustcolor(region_colors, alpha.f = 0.4)  # 非灵长类使用浅色版本

# 自定义颜色映射：灵长类使用深色，非灵长类使用浅色
# color_mapping <- c(
#   rep(colors_primate, each = 1),      # 灵长类颜色
#   rep(colors_nonprimate, each = 1)    # 非灵长类颜色
# )
color_mapping <- c(
  rep(colors_nonprimate, each = 1),     # 非灵长类使用浅色
  rep(colors_primate, each = 1)         # 灵长类使用深色
)
# 绘图
ggplot(data) +
  geom_col(aes(x = Value, y = Region, fill = interaction(Region, Group)),
           width = 0.65, color = "black",
           position = position_dodge(width = 0.70)) +    # 左右分布
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 0参考线
  scale_fill_manual(
    values = color_mapping,
    labels = c("Primate", "NonPrimate")
  ) +
  scale_x_continuous(
    breaks = c(-0.15, 0, 0.3, 0.6, 0.9),  # 手动指定刻度 # 设置横坐标刻度
    limits = c(-0.16, 0.9)                 # 设置横坐标范围
  ) +
  # labs(x = "Feature Value", y = "Brain Region", fill = "Group") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.title.y = element_blank(),  # 去除y轴标题
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# 保存为 PNG 图片
ggsave("pitnum_growth_coupled_SA_expansion_lh.png", width = 3.7, height = 3.7, dpi = 300)


