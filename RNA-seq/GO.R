library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# -----------------------------
# 1. 读取基因列表
# -----------------------------
# 假设你的基因列表在一个文本文件，每行一个基因symbol
gene_file <- "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/RNA-seq/DEGs/LS_down_genes_q0.05_fc0.5_gene_name.txt"
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
write.csv(as.data.frame(ego), "H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/RNA-seq/DEGs/GO/LS_down_GO_enrichment_results.csv", row.names=FALSE)

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
write.csv(as.data.frame(ekegg), "", row.names=FALSE)

# -----------------------------
# 4. 可视化示例
# -----------------------------
# GO气泡图
####LS
ego_sub <- ego
ego_sub@result <- ego@result[c(2, 5, 6, 9, 10, 13, 18, 19, 20, 32, 33, 43, 44, 60, 63, 65, 88, 113, 133), ]

dotplot(ego_sub, showCategory=20) + ggtitle("GO Biological Process Enrichment")


####OS
ego_sub <- ego
ego_sub@result <- ego@result[c(1, 3, 11, 20, 25, 40, 48, 56, 66, 67, 71, 78, 88, 90, 103, 107, 194, 195), ]

dotplot(ego_sub, showCategory=20) + ggtitle("OS GO Biological Process Enrichment")



# KEGG气泡图
dotplot(ekegg, showCategory=20) + ggtitle("KEGG Pathway Enrichment")