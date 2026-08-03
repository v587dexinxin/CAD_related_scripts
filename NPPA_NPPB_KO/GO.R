library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# -----------------------------
# 1. 读取基因列表
# -----------------------------
# 假设你的基因列表在一个文本文件，每行一个基因symbol
gene_file <- "H:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/Confirmation_Experiment/NPPA_NPPB_peak2_敲除/RNA_seq/DEseq2/KO13_up_genes_q0.05_fc1_gene_name.txt"
genes <- readLines(gene_file)
head(genes)

# 将基因symbol转换为Entrez ID
gene_df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_entrez <- gene_df$ENTREZID

# -----------------------------
# 2. GO 富集分析
# -----------------------------
ego <- enrichGO(
  gene          = gene_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",       # 可选: "BP", "MF", "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 输出 GO 富集结果到 CSV
write.csv(as.data.frame(ego), "H:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/Confirmation_Experiment/NPPA_NPPB_peak2_敲除/RNA_seq/DEseq2/WT_VS_NP13_up_GO_enrichment_results.csv", row.names=FALSE)

# -----------------------------
# 3. KEGG 富集分析
# -----------------------------
ekegg <- enrichKEGG(
  gene         = gene_entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05)

# 转换 KEGG ID 为可读 symbol
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# 输出 KEGG 富集结果到 CSV
write.csv(as.data.frame(ekegg), "H:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/Confirmation_Experiment/NPPA_NPPB_peak2_敲除/RNA_seq/DEseq2/WT_VS_NP13_up_KEGG_enrichment_results.csv", row.names=FALSE)

# -----------------------------
# 4. 可视化示例
# -----------------------------
# GO气泡图
dotplot(ego, showCategory=20) + ggtitle("GO Biological Process Enrichment")

# KEGG气泡图
dotplot(ekegg, showCategory=20) + ggtitle("KEGG Pathway Enrichment")