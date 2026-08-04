library(tidyverse)
library(pheatmap)

# 1. 读取 FPKM 矩阵
fpkm <- read.csv(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/RNA-seq/FPKM/HUVEC_FPKM_matrix.csv",
  header = TRUE,
  row.names = 1,
  sep = ",",
  check.names = FALSE
)

# 2. 设置需要展示的基因
genes_show <- c(
  "KLF2", "KLF4", "NOS3", "THBD",
  "NQO1", "HMOX1",
  "ICAM1", "SELE", "CCL2", "CXCL8", "EDN1"
)

# 3. 提取目标基因
genes_found <- intersect(genes_show, rownames(fpkm))
genes_missing <- setdiff(genes_show, rownames(fpkm))

cat("Found genes:\n")
print(genes_found)

cat("Missing genes:\n")
print(genes_missing)

fpkm_sub <- fpkm[genes_found, ]

# 4. log2(FPKM + 1) 转换
fpkm_log2 <- log2(fpkm_sub + 1)

# 5. 可选：按基因做 Z-score，适合展示不同基因在样本间的变化趋势
fpkm_zscore <- t(scale(t(fpkm_log2)))

# 6. 样本分组信息
sample_group <- data.frame(
  Group = c("WT", "WT", "LSS", "LSS", "OSS", "OSS")
)

rownames(sample_group) <- colnames(fpkm_zscore)

# 7. 设置分组颜色
ann_colors <- list(
  Group = c(
    WT  = "#4C78A8",  # blue
    LSS = "#59A14F",  # green
    OSS = "#E15759"   # red
  )
)

# 8. 画热图
pdf("H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/Fig1/Fig1D_selected_genes_FPKM_heatmap.pdf", width = 6, height = 5)
pheatmap(
  fpkm_zscore,
  annotation_col = sample_group,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 10,
  border_color = NA,
  main = "Expression of marker genes"
)


dev.off()


# 9. 保存提取后的表达矩阵
write.table(
  fpkm_sub,
  file = "H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/Fig1/selected_genes_FPKM.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)