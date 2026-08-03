library(DiffBind)
library(DESeq2)


###############LSS_VS_OS
outdir <- "/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/DiffBind/DiffBind_DEseq2"
sample_sheet <- file.path(outdir, "LS_VS_OS_linux.csv")

# -----------------------------
# 1. 读取 DiffBind sample sheet
# -----------------------------
dbObj <- dba(sampleSheet = sample_sheet)

print(dba.show(dbObj))

# -----------------------------
# 2. DiffBind 重新计数
# minOverlap = 1 尽量接近 union peak
# -----------------------------
dbObj <- dba.count(
  dbObj,
  minOverlap = 1,
  bUseSummarizeOverlaps = TRUE
)

# -----------------------------
# 3. 提取 DiffBind count matrix
# -----------------------------
count_gr <- dba.peakset(
  dbObj,
  bRetrieve = TRUE,
  DataType = DBA_DATA_FRAME
)

count_df <- as.data.frame(count_gr)

# 查看列名，确认 peak 坐标列和样本 count 列
print(colnames(count_df))

write.table(
  count_df,
  file = file.path(outdir, "HUVEC_LS_VS_OS_DiffBind_raw_counts_for_DESeq2.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 4. 自动识别 peak 坐标列
# 兼容 seqnames/start/end 或 Chr/Start/End 格式
# -----------------------------
if (all(c("CHR", "START", "END") %in% colnames(count_df))) {
  chr_col <- "CHR"
  start_col <- "START"
  end_col <- "END"
} else if (all(c("Chr", "Start", "End") %in% colnames(count_df))) {
  chr_col <- "Chr"
  start_col <- "Start"
  end_col <- "End"
} else if (all(c("chr", "start", "end") %in% colnames(count_df))) {
  chr_col <- "chr"
  start_col <- "start"
  end_col <- "end"
} else {
  stop("没有找到标准 peak 坐标列，请先检查 colnames(count_df)")
}

# -----------------------------
# 5. 构建 peak annotation
# -----------------------------
peak_anno <- data.frame(
  chr = count_df[[chr_col]],
  start = count_df[[start_col]],
  end = count_df[[end_col]]
)

peak_anno$peak_id <- paste(
  peak_anno$chr,
  peak_anno$start,
  peak_anno$end,
  sep = "_"
)

# -----------------------------
# 6. 提取样本 count 列
# 根据你的样本名包含 LS / OS / LSS / OSS 来识别
# -----------------------------
sample_cols <- grep("LSS|OSS|LS|OS", colnames(count_df), value = TRUE)

print(sample_cols)

if (length(sample_cols) < 4) {
  stop("识别到的样本列太少，请检查 sample_cols 是否正确")
}

count_mat <- count_df[, sample_cols]

# 转为 matrix
count_mat <- as.matrix(count_mat)

# 确保是 integer raw counts
mode(count_mat) <- "integer"

# 设置 peak ID
rownames(count_mat) <- peak_anno$peak_id
rownames(peak_anno) <- peak_anno$peak_id

# -----------------------------
# 7. 构建样本分组信息
# -----------------------------
condition <- ifelse(
  grepl("LSS|LS", sample_cols),
  "LSS",
  "OSS"
)

condition <- factor(condition, levels = c("LSS", "OSS"))

coldata <- data.frame(
  row.names = sample_cols,
  condition = condition
)

print(coldata)

# -----------------------------
# 8. 构建 DESeq2 对象
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

# -----------------------------
# 9. 过滤低 count peaks
# 至少 2 个样本 count >= 10
# -----------------------------
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# 同步过滤 peak annotation
peak_anno_filter <- peak_anno[rownames(dds), ]

# -----------------------------
# 10. DESeq2 差异分析
# -----------------------------
dds <- DESeq(dds)

# OSS vs LSS
res_OSS_vs_LSS <- results(
  dds,
  contrast = c("condition", "OSS", "LSS")
)

# -----------------------------
# 11. 整理结果，并合并 peak 坐标
# -----------------------------
res_df <- as.data.frame(res_OSS_vs_LSS)
res_df$peak_id <- rownames(res_df)

# 合并坐标
res_df <- merge(
  peak_anno_filter,
  res_df,
  by = "peak_id",
  all.y = TRUE
)

# 调整列顺序
res_df <- res_df[, c(
  "peak_id",
  "chr",
  "start",
  "end",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

# -----------------------------
# 12. 添加方向注释
# log2FoldChange > 0 表示 OSS higher
# log2FoldChange < 0 表示 LSS higher
# -----------------------------
res_df$direction <- "not_significant"

res_df$direction[
  !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    res_df$log2FoldChange > 0.5
] <- "OSS_specific"

res_df$direction[
  !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    res_df$log2FoldChange < -0.5
] <- "LSS_specific"

# -----------------------------
# 13. 输出所有结果
# -----------------------------
write.table(
  res_df,
  file = file.path(outdir, "HUVEC_LS_VS_OS_ATAC_deseq2_RLE_BACKGROUND_minOverlap1_all.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 14. 分别输出 LSS / OSS specific peaks
# -----------------------------
LSS_specific <- res_df[
  res_df$direction == "LSS_specific",
]

OSS_specific <- res_df[
  res_df$direction == "OSS_specific",
]

no_sign <- res_df[
  res_df$direction == "not_significant",
]


write.table(
  LSS_specific,
  file = file.path(outdir, "HUVEC_LS_VS_OS_ATAC_DESeq2_LSS_specific_peaks_p0.05_fc0.5.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  OSS_specific,
  file = file.path(outdir, "HUVEC_LS_VS_OS_ATAC_DESeq2_OSS_specific_peaks_p0.05_fc0.5.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  no_sign,
  file = file.path(outdir, "HUVEC_LS_VS_OS_ATAC_DESeq2_common_peaks_p0.05_fc0.5.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)



cat("Total tested peaks:", nrow(res_df), "\n")
cat("LSS specific peaks:", nrow(LSS_specific), "\n")
cat("OSS specific peaks:", nrow(OSS_specific), "\n")
cat("common peaks:", nrow(no_sign), "\n")







































###############WT_VS_OS

library(DiffBind)
library(DESeq2)

outdir <- "/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/DiffBind/DiffBind_DEseq2/WT_VS_OS"
sample_sheet <- file.path(outdir, "WT_VS_OS_linux.csv")

# -----------------------------
# 1. 读取 DiffBind sample sheet
# -----------------------------
dbObj <- dba(sampleSheet = sample_sheet)

print(dba.show(dbObj))

# -----------------------------
# 2. DiffBind 重新计数
# minOverlap = 1 尽量接近 union peak
# -----------------------------
dbObj <- dba.count(
  dbObj,
  minOverlap = 1,
  bUseSummarizeOverlaps = TRUE
)

# -----------------------------
# 3. 提取 DiffBind count matrix
# -----------------------------
count_gr <- dba.peakset(
  dbObj,
  bRetrieve = TRUE,
  DataType = DBA_DATA_FRAME
)

count_df <- as.data.frame(count_gr)

# 查看列名，确认 peak 坐标列和样本 count 列
print(colnames(count_df))

write.table(
  count_df,
  file = file.path(outdir, "HUVEC_WT_VS_OS_DiffBind_raw_counts_for_DESeq2.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 4. 自动识别 peak 坐标列
# 兼容 seqnames/start/end 或 Chr/Start/End 格式
# -----------------------------
if (all(c("CHR", "START", "END") %in% colnames(count_df))) {
  chr_col <- "CHR"
  start_col <- "START"
  end_col <- "END"
} else if (all(c("Chr", "Start", "End") %in% colnames(count_df))) {
  chr_col <- "Chr"
  start_col <- "Start"
  end_col <- "End"
} else if (all(c("chr", "start", "end") %in% colnames(count_df))) {
  chr_col <- "chr"
  start_col <- "start"
  end_col <- "end"
} else {
  stop("没有找到标准 peak 坐标列，请先检查 colnames(count_df)")
}

# -----------------------------
# 5. 构建 peak annotation
# -----------------------------
peak_anno <- data.frame(
  chr = count_df[[chr_col]],
  start = count_df[[start_col]],
  end = count_df[[end_col]]
)

peak_anno$peak_id <- paste(
  peak_anno$chr,
  peak_anno$start,
  peak_anno$end,
  sep = "_"
)

# -----------------------------
# 6. 提取样本 count 列
# 根据你的样本名包含 WT / OS / WT / OSS 来识别
# -----------------------------
sample_cols <- grep("WT|OSS|WT|OS", colnames(count_df), value = TRUE)

print(sample_cols)

if (length(sample_cols) < 4) {
  stop("识别到的样本列太少，请检查 sample_cols 是否正确")
}

count_mat <- count_df[, sample_cols]

# 转为 matrix
count_mat <- as.matrix(count_mat)

# 确保是 integer raw counts
mode(count_mat) <- "integer"

# 设置 peak ID
rownames(count_mat) <- peak_anno$peak_id
rownames(peak_anno) <- peak_anno$peak_id

# -----------------------------
# 7. 构建样本分组信息
# -----------------------------
condition <- ifelse(
  grepl("WT|WT", sample_cols),
  "WT",
  "OSS"
)

condition <- factor(condition, levels = c("WT", "OSS"))

coldata <- data.frame(
  row.names = sample_cols,
  condition = condition
)

print(coldata)

# -----------------------------
# 8. 构建 DESeq2 对象
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

# -----------------------------
# 9. 过滤低 count peaks
# 至少 2 个样本 count >= 10
# -----------------------------
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# 同步过滤 peak annotation
peak_anno_filter <- peak_anno[rownames(dds), ]

# -----------------------------
# 10. DESeq2 差异分析
# -----------------------------
dds <- DESeq(dds)

# OSS vs WT
res_OSS_vs_WT <- results(
  dds,
  contrast = c("condition", "OSS", "WT")
)

# -----------------------------
# 11. 整理结果，并合并 peak 坐标
# -----------------------------
res_df <- as.data.frame(res_OSS_vs_WT)
res_df$peak_id <- rownames(res_df)

# 合并坐标
res_df <- merge(
  peak_anno_filter,
  res_df,
  by = "peak_id",
  all.y = TRUE
)

# 调整列顺序
res_df <- res_df[, c(
  "peak_id",
  "chr",
  "start",
  "end",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

# -----------------------------
# 12. 添加方向注释
# log2FoldChange > 0 表示 OSS higher
# log2FoldChange < 0 表示 WT higher
# -----------------------------
res_df$direction <- "not_significant"

res_df$direction[
  !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    res_df$log2FoldChange > 0.5
] <- "OSS_specific"

res_df$direction[
  !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    res_df$log2FoldChange < -0.5
] <- "WT_specific"

# -----------------------------
# 13. 输出所有结果
# -----------------------------
write.table(
  res_df,
  file = file.path(outdir, "HUVEC_WT_VS_OS_ATAC_deseq2_RLE_BACKGROUND_minOverlap1_all.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 14. 分别输出 WT / OSS specific peaks
# -----------------------------
WT_specific <- res_df[
  res_df$direction == "WT_specific",
]

OSS_specific <- res_df[
  res_df$direction == "OSS_specific",
]

no_sign <- res_df[
  res_df$direction == "not_significant",
]

write.table(
  WT_specific,
  file = file.path(outdir, "HUVEC_WT_VS_OS_ATAC_DESeq2_WT_specific_peaks_p0.05_fc0.5.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  OSS_specific,
  file = file.path(outdir, "HUVEC_WT_VS_OS_ATAC_DESeq2_OSS_specific_peaks_p0.05_fc0.5.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  no_sign,
  file = file.path(outdir, "HUVEC_WT_VS_OS_ATAC_DESeq2_common_peaks_p0.05_fc0.5.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

cat("Total tested peaks:", nrow(res_df), "\n")
cat("WT specific peaks:", nrow(WT_specific), "\n")
cat("OSS specific peaks:", nrow(OSS_specific), "\n")
cat("common peaks:", nrow(no_sign), "\n")









