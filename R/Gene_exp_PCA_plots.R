library(DESeq2)
library(ggplot2)
library(ggrepel)

# 2. 读取 raw count 矩阵
count_data <- read.csv(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/RNA-seq/reads_count/PCA/union_all_reads_count.csv",
  header = TRUE,
  row.names = 1,
  sep = ",",
  check.names = FALSE
)

# 3. 读取样本分组信息
sample_info <- read.table(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/RNA-seq/reads_count/PCA/sample_info.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# 4. 确保 count 矩阵列名与 sample_info 行名一致
count_data <- count_data[, rownames(sample_info)]

# 5. 过滤低表达基因
# 至少在 2 个样本中 count >= 10
keep <- rowSums(count_data >= 10) >= 2
count_data_filtered <- count_data[keep, ]

cat("Genes before filtering:", nrow(count_data), "\n")
cat("Genes after filtering:", nrow(count_data_filtered), "\n")

# 6. 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData = sample_info,
  design = ~ group
)

# 7. vst 标准化
vsd <- vst(dds, blind = TRUE)

# 8. 提取标准化后的表达矩阵
expr_vst <- assay(vsd)

# 9. 做 PCA
pca_res <- prcomp(t(expr_vst), scale. = FALSE)

# 10. 计算 PC1 和 PC2 解释的方差比例
pca_var <- pca_res$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 2)

# 11. 整理 PCA 结果
pca_df <- data.frame(
  sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = sample_info$group
)

# 12. 绘制 PCA 图
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = sample), size = 4) +
  stat_ellipse(
    aes(fill = group),
    geom = "polygon",
    alpha = 0.15,
    color = NA
  ) +
  scale_color_manual(values = c(
    "WT"  = "#1F77B4",  # blue
    "LSS" = "#2CA02C",  # green
    "OSS" = "#D62728"   # red
  )) +
  scale_fill_manual(values = c(
    "WT"  = "#1F77B4",
    "LSS" = "#2CA02C",
    "OSS" = "#D62728"
  )) +
  labs(
    title = "PCA analysis of RNA-seq samples",
    x = paste0("PC1: ", pca_var_percent[1], "% variance"),
    y = paste0("PC2: ", pca_var_percent[2], "% variance"),
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# 13. 显示图形
print(p)

# 14. 保存 PCA 图
ggsave(
  filename = "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/Fig1/Fig1C_RNAseq_PCA.pdf",
  plot = p,
  width = 6,
  height = 5
)



# 15. 保存 PCA 坐标结果
write.table(
  pca_df,
  file = "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/Fig1/Fig1C_RNAseq_PCA_coordinates.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)