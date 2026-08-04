library(DESeq2)

# 读取 count 矩阵
count_df <- read.table(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/STARR-seq/DiffBind/Regulatory_Element/edgR/fc_0/DAVs/HUVEC_STARR_seq_cDNA.txt",
  header = FALSE,
  sep = "\t",
  check.names = FALSE,
  skip = 2,
  col.names = c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'HUVEC_R1_nonrisk', 'HUVEC_R2_nonrisk', 'HUVEC_R1_risk', 'HUVEC_R2_risk')
)

# 查看前几行
head(count_df)

# 提取注释信息
anno_df <- count_df[, c("Geneid", "Chr", "Start", "End", "Strand", "Length")]

# 提取 count 数据
count_mat <- count_df[, c(
  "HUVEC_R1_nonrisk",
  "HUVEC_R2_nonrisk",
  "HUVEC_R1_risk",
  "HUVEC_R2_risk"
)]

# 设置 rownames
rownames(count_mat) <- count_df$Geneid

# 确保 count 是整数
count_mat <- round(as.matrix(count_mat))
mode(count_mat) <- "integer"

# 构建样本信息
coldata <- data.frame(
  row.names = colnames(count_mat),
  condition = c("nonrisk", "nonrisk", "risk", "risk")
)

# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

# 过滤低 count 区域，可根据需要调整
dds <- dds[rowSums(counts(dds)) > 10, ]

# DESeq2 归一化
dds <- estimateSizeFactors(dds)

# 提取归一化后的 count
norm_counts <- counts(dds, normalized = TRUE)

# 查看 size factors
sizeFactors(dds)


# 转成 data.frame
norm_counts_df <- as.data.frame(norm_counts)

# 添加 Geneid
norm_counts_df$Geneid <- rownames(norm_counts_df)

# 合并注释信息
result_df <- merge(
  anno_df,
  norm_counts_df,
  by = "Geneid"
)

# 按原始顺序整理
result_df <- result_df[, c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length",
  "HUVEC_R1_nonrisk",
  "HUVEC_R2_nonrisk",
  "HUVEC_R1_risk",
  "HUVEC_R2_risk"
)]

# 保存结果
write.table(
  result_df,
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/STARR-seq/DiffBind/Regulatory_Element/edgR/fc_0/DAVs/HUVEC_STARR_seq_cDNA_DESeq2_normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





#############################plasmid
# 读取 count 矩阵
count_df <- read.table(
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/STARR-seq/DiffBind/Regulatory_Element/edgR/fc_0/DAVs/HUVEC_STARR_seq_plasmid.txt",
  header = FALSE,
  sep = "\t",
  check.names = FALSE,
  skip = 2,
  col.names = c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'HUVEC_R1_nonrisk', 'HUVEC_R2_nonrisk', 'HUVEC_R1_risk', 'HUVEC_R2_risk')
)

# 查看前几行
head(count_df)

# 提取注释信息
anno_df <- count_df[, c("Geneid", "Chr", "Start", "End", "Strand", "Length")]

# 提取 count 数据
count_mat <- count_df[, c(
  "HUVEC_R1_nonrisk",
  "HUVEC_R2_nonrisk",
  "HUVEC_R1_risk",
  "HUVEC_R2_risk"
)]

# 设置 rownames
rownames(count_mat) <- count_df$Geneid

# 确保 count 是整数
count_mat <- round(as.matrix(count_mat))
mode(count_mat) <- "integer"

# 构建样本信息
coldata <- data.frame(
  row.names = colnames(count_mat),
  condition = c("nonrisk", "nonrisk", "risk", "risk")
)

# 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

# 过滤低 count 区域，可根据需要调整
dds <- dds[rowSums(counts(dds)) > 10, ]

# DESeq2 归一化
dds <- estimateSizeFactors(dds)

# 提取归一化后的 count
norm_counts <- counts(dds, normalized = TRUE)

# 查看 size factors
sizeFactors(dds)


# 转成 data.frame
norm_counts_df <- as.data.frame(norm_counts)

# 添加 Geneid
norm_counts_df$Geneid <- rownames(norm_counts_df)

# 合并注释信息
result_df <- merge(
  anno_df,
  norm_counts_df,
  by = "Geneid"
)

# 按原始顺序整理
result_df <- result_df[, c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length",
  "HUVEC_R1_nonrisk",
  "HUVEC_R2_nonrisk",
  "HUVEC_R1_risk",
  "HUVEC_R2_risk"
)]

# 保存结果
write.table(
  result_df,
  "H:/work/Postdoctoral/GWAS疾病位点检测/results/STARR-seq/DiffBind/Regulatory_Element/edgR/fc_0/DAVs/HUVEC_STARR_seq_plasmid_DESeq2_normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)








