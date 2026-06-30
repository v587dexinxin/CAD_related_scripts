library(DESeq2)
Data <- read.table("H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/reads_count/Union_Reads_count_norm.txt", header=T, row.names=1, sep="\t")
sample <- read.table("H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/DESeq2/HCT116_cDNA_VS_plasmid.csv", header=T, row.names=1, com='', quote='', check.names=F, sep=",", colClasses="factor")
data <- Data[c('HCT116.R1.4_S1','HCT116.R2.8_S2','HCT116.P.R1.1_S5','HCT116.P.R2.2_S6')]
data <- data[rowSums(data)>2,]

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                            colData = sample,  design= ~ conditions)

dds <- DESeq(ddsFullCountTable)

rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)


sampleA <- 'cDNA'
sampleB <- 'plasmid'


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


write.table(res, file='H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/DESeq2/HCT116_deseq2_norm.csv', sep=",", quote=F, row.names=F)


####HUVEC
Data <- read.table("H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/reads_count/Union_Reads_count_norm.txt", header=T, row.names=1, sep="\t")
sample <- read.table("H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/DESeq2/HUVEC_cDNA_VS_plasmid_samplesheet.csv", header=T, row.names=1, com='', quote='', check.names=F, sep=",", colClasses="factor")
data <- Data[c('HUVEC.R1.3_S3','HUVEC.R2.6_S4','HUVEC.P.R1.5_S7','HUVEC.P.R2.7_S8')]
data <- data[rowSums(data)>2,]

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                            colData = sample,  design= ~ conditions)

dds <- DESeq(ddsFullCountTable)

rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)


sampleA <- 'cDNA'
sampleB <- 'plasmid'


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


write.table(res, file='H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/DESeq2/HUVEC_deseq2_norm.csv', sep=",", quote=F, row.names=F)








