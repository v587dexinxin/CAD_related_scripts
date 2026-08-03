library(DESeq2)

count_df <- read.table(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/one-dimensional_new/bam/HiCoatis_LS_OS_union_peaks_counts.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

# featureCounts 前几列是注释信息，后面是样本 counts
counts <- count_df[, 7:ncol(count_df)]
rownames(counts) <- count_df$Geneid
colnames(counts) <- c("LSS_R1", "LSS_R2", "LSS_R3", "OSS_R1" , "OSS_R2" , "OSS_R3")

# 样本信息
condition <- factor(c("LSS", "LSS", "LSS", "OSS", "OSS", "OSS"))
coldata <- data.frame(row.names = colnames(counts), condition = condition)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# 去除低 counts peaks
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# 查看 size factors
sizeFactors(dds)

# 提取 peak 位置信息
peak_anno <- count_df[, c("Geneid", "Chr", "Start", "End", "Strand", "Length")]

# 获得归一化后的 counts
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$Geneid <- rownames(norm_counts_df)

# normalized counts 加入 peak 信息
norm_counts_anno <- merge(
  peak_anno,
  norm_counts_df,
  by = "Geneid",
  all.y = TRUE
)

write.table(
  norm_counts_anno,
  file = "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/one-dimensional_new/bam/HiCoatis_LS_OS_union_peaks_DESeq2_normalized_counts_with_peak_info.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 差异分析
res <- results(dds, contrast = c("condition", "LSS", "OSS"))

res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)

# DESeq2 结果加入 peak 信息
res_df_anno <- merge(
  peak_anno,
  res_df,
  by = "Geneid",
  all.y = TRUE
)

# 按 padj 排序
res_df_anno <- res_df_anno[order(res_df_anno$padj), ]

write.table(
  res_df_anno,
  file = "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/one-dimensional_new/bam/HiCoatis_LSS_vs_OSS_DESeq2_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


sizeFactors(dds)

boxplot(
  log2(counts + 1),
  las = 2,
  outline = FALSE,
  ylab = "log2(counts + 1)"
)

colSums(counts)



#########HiCoatis_union_peaks_MA_plot


library(ggplot2)

# counts 是 union peaks 的 raw count matrix
# 样本列名请根据你的实际列名调整LSS_R1OSS_R1
lss_cols <- c("LSS_R1",
              "LSS_R2",
              "LSS_R3")

oss_cols <- c("OSS_R1",
              "OSS_R2",
              "OSS_R3")

df <- data.frame(
  LSS_mean = rowMeans(counts[, lss_cols]),
  OSS_mean = rowMeans(counts[, oss_cols])
)

df$A <- log2((df$LSS_mean + df$OSS_mean) / 2 + 1)
df$M <- log2((df$OSS_mean + 1) / (df$LSS_mean + 1))

p <- ggplot(df, aes(x = A, y = M)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Average signal: log2(mean counts + 1)",
    y = "log2(OSS / LSS)",
    title = "MA plot of Hi-Coatis signal over union peaks"
  )

print(p)

ggsave("H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_MA_plot.pdf", p, width = 5, height = 4)





library(ggplot2)

res_df$direction <- "Not significant"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "OSS-up"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "LSS-up"

# 需要确保 rownames 和 counts 行名一致
res_df$LSS_mean <- rowMeans(counts[res_df$peak_id, lss_cols])
res_df$OSS_mean <- rowMeans(counts[res_df$peak_id, oss_cols])

res_df$A <- log2((res_df$LSS_mean + res_df$OSS_mean) / 2 + 1)
res_df$M <- log2((res_df$OSS_mean + 1) / (res_df$LSS_mean + 1))

p <- ggplot(res_df, aes(x = A, y = M, color = direction)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "OSS-up" = "#D73027",
      "LSS-up" = "#4575B4",
      "Not significant" = "grey70"
    )
  ) +
  theme_classic() +
  labs(
    x = "Average signal: log2(mean counts + 1)",
    y = "log2(OSS / LSS)",
    color = NULL,
    title = "Differential Hi-Coatis signal over union peaks"
  )

print(p)

ggsave("H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_MA_plot_DESeq2_colored.pdf", p, width = 5.5, height = 5)






#######LS_OS信号强度分布
df$LSS_log2 <- log2(df$LSS_mean + 1)
df$OSS_log2 <- log2(df$OSS_mean + 1)

plot_df <- data.frame(
  Signal = c(df$LSS_log2, df$OSS_log2),
  Group = rep(c("LSS", "OSS"), each = nrow(df))
)

p2 <- ggplot(plot_df, aes(x = Signal, fill = Group)) +
  geom_density(alpha = 0.35) +
  theme_classic() +
  labs(
    x = "log2(mean counts + 1)",
    y = "Density",
    title = "Distribution of Hi-Coatis signal over union peaks"
  )

print(p2)

ggsave("H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/S1/HiCoatis_union_peaks_density.pdf", p2, width = 5, height = 4)


















