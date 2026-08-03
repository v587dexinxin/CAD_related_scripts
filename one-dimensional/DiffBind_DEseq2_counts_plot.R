library(DESeq2)
library(ggplot2)
library(reshape2)
library(dplyr)

# 读取 featureCounts 输出
count_df <- read.table(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/one-dimensional_new/bam/HiCoatis_union_all_peaks_counts.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

# 提取 count 矩阵
# featureCounts 前 6 列通常是 Geneid, Chr, Start, End, Strand, Length
counts <- count_df[, 7:ncol(count_df)]
rownames(counts) <- count_df$Geneid

# 修改样本名，根据你的截图
colnames(counts) <- c(
  "WT_Rep1", "WT_Rep2", "WT_Rep3",
  "LSS_Rep1", "LSS_Rep2", "LSS_Rep3",
  "OSS_Rep1", "OSS_Rep2", "OSS_Rep3"
)

# 样本分组信息
condition <- factor(
  c(rep("WT", 3), rep("LSS", 3), rep("OSS", 3)),
  levels = c("WT", "LSS", "OSS")
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)




#############DESeq2_归一化
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = coldata,
  design = ~ condition
)

# 可选：过滤极低 count peaks
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 估计 size factors
dds <- estimateSizeFactors(dds)

# 提取归一化后的 counts
norm_counts <- counts(dds, normalized = TRUE)

# log2 转换
log2_norm_counts <- log2(norm_counts + 1)



sizeFactors(dds)

###############归一化后 log2(counts + 1) 箱线图
plot_df <- as.data.frame(log2_norm_counts)
plot_df$Peak <- rownames(plot_df)

# 提取样本矩阵
box_mat <- plot_df[, setdiff(colnames(plot_df), "Peak")]

# 固定样本顺序
sample_order <- c(
  "WT_Rep1", "WT_Rep2", "WT_Rep3",
  "LSS_Rep1", "LSS_Rep2", "LSS_Rep3",
  "OSS_Rep1", "OSS_Rep2", "OSS_Rep3"
)

box_mat <- box_mat[, sample_order]

# 设置颜色
box_cols <- c(
  rep("#4C78A8", 3),
  rep("#54A24B", 3),
  rep("#E45756", 3)
)


pdf(
  "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_log2_normalized_counts_boxplot_baseR_outline.pdf",
  width = 7,
  height = 5
)

boxplot(
  box_mat,
  col = "white",
  border = box_cols,
  outline = FALSE,
  boxwex = 0.6,
  las = 2,
  main = "Normalized Hi-Coatis signal across union peaks",
  ylab = "log2(normalized counts + 1)",
  xlab = "",
  cex.axis = 0.9,
  cex.lab = 1.1,
  cex.main = 1.2
)

legend(
  "right",
  legend = c("WT", "LSS", "OSS"),
  col = c("#4C78A8", "#54A24B", "#E45756"),
  lwd = 2,
  bty = "n"
)

dev.off()







###############HiCoatis_union_peaks_signal_density_WT_LSS_OSS
# 计算每组 mean normalized counts，再 log2(mean + 1)
signal_mean_df <- data.frame(
  Peak = rownames(norm_counts),
  WT = log2(rowMeans(norm_counts[, c("WT_Rep1", "WT_Rep2", "WT_Rep3")], na.rm = TRUE) + 1),
  LSS = log2(rowMeans(norm_counts[, c("LSS_Rep1", "LSS_Rep2", "LSS_Rep3")], na.rm = TRUE) + 1),
  OSS = log2(rowMeans(norm_counts[, c("OSS_Rep1", "OSS_Rep2", "OSS_Rep3")], na.rm = TRUE) + 1)
)

# 转换为长格式
signal_long <- melt(
  signal_mean_df,
  id.vars = "Peak",
  variable.name = "Condition",
  value.name = "log2_mean_norm_count"
)

signal_long$Condition <- factor(
  signal_long$Condition,
  levels = c("WT", "LSS", "OSS")
)

# 绘制 density plot
p_density <- ggplot(
  signal_long,
  aes(x = log2_mean_norm_count, color = Condition, fill = Condition)
) +
  geom_density(
    alpha = 0.25,
    linewidth = 1
  ) +
  scale_color_manual(values = c(
    "WT" = "#4C78A8",
    "LSS" = "#54A24B",
    "OSS" = "#E45756"
  )) +
  scale_fill_manual(values = c(
    "WT" = "#4C78A8",
    "LSS" = "#54A24B",
    "OSS" = "#E45756"
  )) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Hi-Coatis signal density across union peaks",
    x = "log2(mean normalized counts + 1)",
    y = "Density"
  )

p_density

ggsave(
  "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_signal_density_WT_LSS_OSS.pdf",
  p_density,
  width = 6,
  height = 5
)








############HiCoatis_union_peaks_MA_plot

ma_df <- data.frame(
  Peak = rownames(log2_norm_counts),
  
  WT_mean  = rowMeans(log2_norm_counts[, c("WT_Rep1", "WT_Rep2", "WT_Rep3")], na.rm = TRUE),
  LSS_mean = rowMeans(log2_norm_counts[, c("LSS_Rep1", "LSS_Rep2", "LSS_Rep3")], na.rm = TRUE),
  OSS_mean = rowMeans(log2_norm_counts[, c("OSS_Rep1", "OSS_Rep2", "OSS_Rep3")], na.rm = TRUE)
)

# OSS vs LSS
ma_df <- data.frame(
  Peak = rownames(log2_norm_counts),
  LSS_mean = rowMeans(log2_norm_counts[, c("LSS_Rep1", "LSS_Rep2", "LSS_Rep3")], na.rm = TRUE),
  OSS_mean = rowMeans(log2_norm_counts[, c("OSS_Rep1", "OSS_Rep2", "OSS_Rep3")], na.rm = TRUE)
)

# MA plot 数据
# A: 平均信号强度
# M: OSS 相对于 LSS 的变化
ma_df <- ma_df %>%
  mutate(
    A = (LSS_mean + OSS_mean) / 2,
    M = LSS_mean - OSS_mean
  )

# 画整体偏移效果，不按颜色区分
p_ma <- ggplot(ma_df, aes(x = A, y = M)) +
  geom_point(
    color = "#B39DDB",
    alpha = 0.45,
    size = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none"
  ) +
  labs(
    title = "Differential Hi-Coatis signal over union peaks",
    x = "Average signal: log2(mean counts + 1)",
    y = "log2(LSS / OSS)"
  )

p_ma


ggsave(
  "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_MA_plot_LSS_vs_OSS_overall_shift_smooth.pdf",
  p_ma,
  width = 5.5,
  height = 5
)













# WT vs LSS
ma_df <- data.frame(
  Peak = rownames(log2_norm_counts),
  LSS_mean = rowMeans(log2_norm_counts[, c("LSS_Rep1", "LSS_Rep2", "LSS_Rep3")], na.rm = TRUE),
  WT_mean = rowMeans(log2_norm_counts[, c("WT_Rep1", "WT_Rep2", "WT_Rep3")], na.rm = TRUE)
)

# MA plot 数据
# A: 平均信号强度
# M: WT 相对于 LSS 的变化
ma_df <- ma_df %>%
  mutate(
    A = (LSS_mean + WT_mean) / 2,
    M = WT_mean - LSS_mean
  )

# 画整体偏移效果，不按颜色区分
p_ma <- ggplot(ma_df, aes(x = A, y = M)) +
  geom_point(
    color = "#B39DDB",
    alpha = 0.45,
    size = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none"
  ) +
  labs(
    title = "Differential Hi-Coatis signal over union peaks",
    x = "Average signal: log2(mean counts + 1)",
    y = "log2(WT / LSS)"
  )

p_ma

ggsave(
  "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_MA_plot_WT_vs_LSS_overall_shift_smooth.pdf",
  p_ma,
  width = 5.5,
  height = 5
)





# WT vs OSS
ma_df <- data.frame(
  Peak = rownames(log2_norm_counts),
  OSS_mean = rowMeans(log2_norm_counts[, c("OSS_Rep1", "OSS_Rep2", "OSS_Rep3")], na.rm = TRUE),
  WT_mean = rowMeans(log2_norm_counts[, c("WT_Rep1", "WT_Rep2", "WT_Rep3")], na.rm = TRUE)
)

# MA plot 数据
# A: 平均信号强度
# M: WT 相对于 OSS 的变化
ma_df <- ma_df %>%
  mutate(
    A = (OSS_mean + WT_mean) / 2,
    M = WT_mean - OSS_mean
  )

# 画整体偏移效果，不按颜色区分
p_ma <- ggplot(ma_df, aes(x = A, y = M)) +
  geom_point(
    color = "#B39DDB",
    alpha = 0.45,
    size = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none"
  ) +
  labs(
    title = "Differential Hi-Coatis signal over union peaks",
    x = "Average signal: log2(mean counts + 1)",
    y = "log2(WT / OSS)"
  )

p_ma



ggsave(
  "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_MA_plot_WT_vs_OSS_overall_shift_smooth.pdf",
  p_ma,
  width = 5.5,
  height = 5
)





