# 清空环境
rm(list = ls())
setwd("Y:/results/Evolution_cortical_shape/statistic/kinship/")

library(R.matlab)
library(vegan)
library(ggplot2)
library(viridis)
library(pheatmap)
library(ecodist)
# 读取数据
DIST_mat <- readMat("DIST_normalized.mat")$DIST.normalized  # 90x90 亲缘关系矩阵，已归一化到 [-1, 1]
SIM_mat <- readMat("pit_corr_Nring_order.mat")$pit.corr.order     # 90x90 pit 相似性矩阵 按照圈数定
# 将相似性矩阵转为距离矩阵：距离 = 1 - 相似性
dist_phylo <- as.dist(1 - DIST_mat)
dist_pit <- as.dist(1 - SIM_mat)

feat <- read.csv("Y:/results/Evolution_cortical_shape/cortical_features/cortical_feature.csv", header = TRUE, stringsAsFactors = FALSE)
# 取表面积与脑体积向量
surf_area <- as.numeric(feat[[2]])
brain_vol <- as.numeric(feat[[3]])

# 用 dist() 直接生成“差值距离”（对一维向量等同于绝对差）
dist_area <- dist(surf_area, method = "euclidean")  # 皮层表面积差
dist_vol  <- dist(brain_vol,  method = "euclidean") # 脑体积差

# 可选：对协变量距离做标准化（让回归系数更可比；不改变显著性本质）
standardize_dist <- function(d) {
  v <- as.vector(d)
  v <- scale(v)                     # z-score
  attr(v, "Size")  <- attr(d, "Size")
  attr(v, "Diag")  <- FALSE
  attr(v, "Upper") <- FALSE
  class(v) <- "dist"
  return(v)
}
dist_area_z <- standardize_dist(dist_area)
dist_vol_z  <- standardize_dist(dist_vol)

# 计算 Mantel 检验
mantel_result <- mantel(dist_pit, dist_phylo, method = "pearson", permutations = 999)
print(mantel_result)

#-----------------------------
# 4) MRM：pit不相似度 ~ 系统发育距离 + 表面积差 + 脑体积差
#    - nperm 建议 >= 9999
#    - 默认置换方案为 QAP 类型（行列同时置换）更稳健
#-----------------------------
set.seed(12345)
mrm_fit <- MRM(dist_pit ~ dist_phylo + dist_area_z + dist_vol_z, nperm = 10000)
print(mrm_fit)

# 结果整理为数据框，便于可视化/汇报
# 从 mrm_fit$coef 构造整洁数据框
coef_mat <- mrm_fit$coef
coef_df <- data.frame(
  term   = rownames(coef_mat),
  coef   = coef_mat[, "dist_pit"],  # 第一列就是系数
  p_perm = coef_mat[, "pval"],
  row.names = NULL
)
coef_df

#-----------------------------
# 5) 简要可视化：系数与置换p值
#-----------------------------
library(ggplot2)
p <- ggplot(subset(coef_df, term != "Int"),
            aes(x = reorder(term, coef), y = coef, fill = p_perm < 0.05)) +
  geom_col(width = 0.65) +
  coord_flip() +
  scale_fill_manual(
    values = c("grey70", "steelblue"),
    labels = c("Not significant", "p_perm < 0.05"),   # 英文图例
    name   = "Significance"                           # 图例标题
  ) +
  labs(
    x = NULL,
    y = "MRM regression coefficient",
    title = "Multiple regression on distance matrices"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",          # 图例放在右边
    legend.title = element_text(size=12),
    legend.text  = element_text(size=10)
  )
ggsave("Y:/results/Evolution_cortical_shape/statistic/kinship/MRM_coefficients.png", p, width = 6, height = 2, dpi = 300)

#-----------------------------
# 6) 可选稳健性：把协变量不标准化再跑一次，比较系数方向与显著性是否一致
#-----------------------------
set.seed(12345)
mrm_fit_raw <- MRM(dist_pit ~ dist_phylo + dist_area + dist_vol, nperm = 10000)
print(mrm_fit_raw)


############ pic2
# 建立映射表
label_map <- c(
  "dist_phylo"  = "Phylo",
  "dist_area_z" = "CSA",
  "dist_vol_z"  = "BV"
)

# 给 coef_df 添加新列 label
coef_df$label <- label_map[coef_df$term]
# 固定顺序
coef_df$label <- factor(coef_df$label, levels = c("Phylo", "CSA", "BV"))
# 画图
p <- ggplot(subset(coef_df, term != "Int"),
            aes(x = label, y = coef, fill = p_perm < 0.05)) +
  geom_col(width = 0.65) +
  scale_fill_manual(
    values = c("grey70", "steelblue"),
    labels = c("Not significant", "p_perm < 0.05"),
    name   = "Significance"
  ) +
  labs(
    x = NULL,
    y = expression(beta),
    title = "MRM"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.95, 0.95),      # 图例放到右上角
    legend.justification = c("right","top"),  # 对齐方式
    legend.background = element_rect(fill = alpha("white", 0.6), colour = NA), # 半透明背景
    legend.title = element_text(size=12),
    legend.text  = element_text(size=10)
  )

ggsave("Y:/results/Evolution_cortical_shape/statistic/kinship/MRM_coefficients.png",
       p, width = 3.3 , height = 3.8, dpi = 300)
