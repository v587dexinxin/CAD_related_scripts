####LS_VS_OS
library(DiffBind)
dbObj <- dba(sampleSheet="/scratch/2025-05-19/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/DiffBind/LS_VS_OS_linux.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("/scratch/2025-05-19/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/DiffBind/LS_VS_OS_PCA.pdf", pointsize=10)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

#dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="/scratch/2025-05-19/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/DiffBind/LS_vs_OS_deseq2.csv", sep=",", quote=F,row.names = FALSE)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="/scratch/2025-05-19/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/DiffBind/LS_vs_OS_new_edgeR.csv", sep=",", quote=F, row.names = FALSE)
