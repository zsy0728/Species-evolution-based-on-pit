# Differences in the rate of evolution of different brain networks
# 
#------------------------------------------------------
rm(list = ls())
dir <- "Y:/sda/results/Evolution_cortical_shape/statistic/evo_speedy/"
setwd(dir)

library(ggplot2)
library(reshape2)
library(RColorBrewer)

library(fmsb)
library("R.matlab")

non <- readMat(file.path("nonpri_yeo7.mat"))
pri <- readMat(file.path("pri_yeo7.mat"))

M_count <- table(pri)
N_count <- table(non)

# 确保编号范围一致（1-7）
brain_labels <- as.character(1:7)
M_freq <- as.numeric(M_count[brain_labels])
N_freq <- as.numeric(N_count[brain_labels])

# 如果有空值，将其设为0
M_freq[is.na(M_freq)] <- 0
N_freq[is.na(N_freq)] <- 0

# 设置脑区名称
region_names <- c("Visual", "SMN", "DAN", "VAN", "Limbic", "FPN", "DMN")

# 创建数据框（用于雷达图）
data_primate <- data.frame(
  Max = rep(max(c(M_freq, N_freq)), length(region_names)),  # 最大值
  Min = rep(0, length(region_names)),                      # 最小值
  Primate = M_freq                                         # 灵长类数据
)
rownames(data_primate) <- region_names
data_primate_t <- as.data.frame(t(data_primate))

data_nonprimate <- data.frame(
  Max = rep(max(c(M_freq, N_freq)), length(region_names)),  # 最大值
  Min = rep(0, length(region_names)),                      # 最小值
  NonPrimate = N_freq                                      # 非灵长类数据
)
rownames(data_nonprimate) <- region_names
data_nonprimate_t <- as.data.frame(t(data_nonprimate))

pdf("YEO7_Evolution_Radar.pdf", width = 10, height = 10)
par(mfrow = c(1, 2)) 

radarchart(
  data_primate_t,
  axistype = 1,
  cglcol = "darkgray",  # 刻度线颜色
  cglty = 1,            # 刻度线类型
  axislabcol = "darkgray",  # 轴标签颜色
  pcol = "#E64B35FF",       # 线条颜色
  plwd = 2,               # 线条宽度
  pfcol = rgb(0.9, 0.3, 0.2, 0.5),  # 填充颜色
  vlcex = 1.2            # 加粗并增大脑区名称字体大小
)
# title("Primate", cex.main = 1.5)

# 绘制非灵长类雷达图
radarchart(
  data_nonprimate_t,
  axistype = 1,
  cglcol = "darkgray",  # 刻度线颜色
  cglty = 1,            # 刻度线类型
  axislabcol = "darkgray",  # 轴标签颜色
  pcol = "#4DBBD5FF",       # 线条颜色
  plwd = 2,               # 线条宽度
  pfcol = rgb(0.3, 0.73, 0.83, 0.5),   # 填充颜色
  vlcex = 1.2
)
# title("Non-Primate", cex.main = 1.5)

# # 添加图例到下方
# par(fig = c(0, 1, 0, 0.2), new = TRUE)  # 设置图例区域
# plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
# legend(
#   "center",  # 图例位置
#   legend = c("Primate", "Non-Primate"),
#   col = c("#E64B35FF", "#4DBBD5FF"),
#   pch = 15,  # 方块
#   pt.cex = 2,  # 图例点大小
#   bty = "n",   # 无边框
#   cex = 1.2    # 图例文本大小
# )

# 关闭 PDF 设备
dev.off()
