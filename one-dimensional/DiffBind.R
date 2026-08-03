####LS_VS_OS
library(DiffBind)
dbObj <- dba(sampleSheet="/scratch/2026-05-11/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/DiffBind/new/LS_VS_OS_linux_1.csv")
dbObj <- dba.count(dbObj, minOverlap=1, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)


##########background_norm
dbObj <- dba.contrast(
  dbObj,
  categories = DBA_FACTOR,
  minMembers = 2
)

dbObj <- dba.normalize(
  dbObj,
  method = DBA_DESEQ2,
  normalize = DBA_NORM_RLE,
  library = DBA_LIBSIZE_FULL,
  background = TRUE
)

dbObj <- dba.analyze(
  dbObj,
  method = DBA_DESEQ2
)

pdf("HUVEC_LS_VS_OS_MA_DESeq2_RLE_FULL_BACKGROUND.pdf")
dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = 1)
dev.off()

res <- dba.report(
  dbObj,
  method = DBA_DESEQ2,
  contrast = 1,
  th = 1
)
###############

pdf("/scratch/2025-12-15/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/DiffBind/HUVEC_LS_VS_OS_PCA.pdf", pointsize=10)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

#dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="/scratch/2025-12-15/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/DiffBind/HUVEC_LS_VS_OS_deseq2.csv", sep=",", quote=F,row.names = FALSE)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="/scratch/2025-12-15/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/DiffBind/HUVEC_LS_VS_OS_edgeR.csv", sep=",", quote=F, row.names = FALSE)
