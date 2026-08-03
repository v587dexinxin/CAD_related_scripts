library(DESeq2)

#########WT_VS_NP13
Data <- read.table("H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/STARR-seq/verification_experiments/Confirmation_Experiment/NPPA_NPPB_peak2_ÇÃ³ý/RNA_seq/reads_count/NPPA_NPPB_peak2_KO_union_all_reads_count.csv", header=T, row.names=1, sep=",")
sample <- read.table("H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/STARR-seq/verification_experiments/Confirmation_Experiment/NPPA_NPPB_peak2_ÇÃ³ý/RNA_seq/reads_count/WT_VS_NP13.csv", header=T, row.names=1, com='', quote='', check.names=F, sep=",", colClasses="factor")
data <- Data[c('WT_116_R1_Count' , 'WT_116_R2_Count' , 'NP_13_R1_Count' , 'NP_13_R2_Count')]
data <- data[rowSums(data)>2,]

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                            colData = sample,  design= ~ conditions)

dds <- DESeq(ddsFullCountTable)

rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)


sampleA <- 'WT'
sampleB <- 'NP13'


contrastV <- c("conditions", sampleA, sampleB)
res <- results(dds,  contrast=contrastV)

baseA <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleA]

if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
head(baseMeanA)


baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
head(baseMeanB)


res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
res <- cbind(Gene_Name=rownames(res), as.data.frame(res))
res$padj[is.na(res$padj)] <- 1

head(res)
res <- res[order(res$padj),]


write.table(res, file='H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/STARR-seq/verification_experiments/Confirmation_Experiment/NPPA_NPPB_peak2_ÇÃ³ý/RNA_seq/DEseq2/HCT116_WT_VS_NP13_KO_deseq2.csv', sep=",", quote=F, row.names=F)
