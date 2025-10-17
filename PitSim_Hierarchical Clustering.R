# 可视化没有做层次聚类的时候五个label的色块
# 加载必要包
library(R.matlab)
library(readxl)
library(tidyverse)
library(cowplot)

# 加载数据
order_idx <- readMat("phytree_indices.mat")$phytree.indices |> as.vector()

infodir <- 'Y:/sda/results/Evolution_cortical_shape/data_info/'
name_list <- read.table(file.path(infodir, "phyTree_name_list.txt"), header = FALSE, stringsAsFactors = FALSE)
node_ages <- read_excel(file.path(infodir, "node_ages.xlsx"))
species_info <- read_excel(file.path(infodir, "Species_info.xlsx"))

# 按照 phylogeny 排序
species_info <- species_info[order_idx, ]

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
ggsave("species_blocks_plot.pdf", p_main, width = 20, height = 4)


# color code
library(ggplot2)

# 按指定顺序定义颜色（名称仅用于内部区分，不显示）
color_order <- list(
  Primata = "#E64B35FF",
  Other = "#4DBBD5FF",
  Arboreal = "#b39cd0",
  Arboreal_alpha = adjustcolor("#b39cd0", alpha.f = 0.5),
  Terrestrial = "#6b9b72",
  Terrestrial_alpha = adjustcolor("#6b9b72", alpha.f = 0.5),
  Fossorial = "#2d8ba8",
  Fossorial_alpha = adjustcolor("#2d8ba8", alpha.f = 0.5),
  Social = "#f0a0e9",
  Social_alpha = adjustcolor("#f0a0e9", alpha.f = 0.5)
)
# 转换为数据框，添加位置索引（确保纵向排布）
color_df <- data.frame(
  color = unlist(color_order),
  y_pos = 1:length(color_order)  # y轴位置从1到10，实现纵向排列
)
# 绘制纵向排布的正方形色块
ggplot(color_df, aes(x = 1, y = y_pos, fill = color)) +
  geom_tile(width = 1, height = 1, color = "white", size = 1) +  # 正方形色块，白色边框
  scale_fill_identity() +  # 使用实际颜色
  scale_y_reverse(limits = c(max(color_df$y_pos) + 0.5, 0.5)) +  # 反转y轴，使第一个元素在顶部
  theme_void() +  # 去除所有边框、坐标轴和文字
  theme(
    plot.margin = margin(20, 20, 20, 20, "pt"),  # 整体边距
    panel.background = element_rect(fill = "white", color = NA)  # 白色背景
  ) +
  coord_fixed(ratio = 1)  # 强制正方形比例，确保色块为正方块
# 保存图形（确保正方形比例）
ggsave(
  "Y:/results/Evolution_cortical_shape/figure/color_squares.png",
  width = 3,
  height = 12,
  dpi = 300,
  bg = "white"  # 白色背景
)

