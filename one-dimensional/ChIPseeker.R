library(ChIPseeker)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("org.Mm.eg.db")
library("GenomicFeatures")

spompe <- makeTxDbFromGFF('H:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')
files <- list(WT = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/peaks/union_peaks/HiRPC_WT_allreps_q0.05_peaks_sorted_merged.bed') , 
              LS = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/peaks/union_peaks/HiRPC_LS_allreps_q0.05_peaks_sorted_merged.bed') , 
              OS = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/peaks/union_peaks/HiRPC_OS_allreps_q0.05_peaks_sorted_merged.bed'))


##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 2000) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)


pdf('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/peaks/union_peaks/HUVEC_WT_VS_LS_OS_union_q0.05_peaks_anno.pdf')
plotAnnoBar(peakAnnoList)
dev.off()







###########specific_peaks

files <- list(LS_speci = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/HUVEC_LS_VS_OS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed') , 
              OS_speci = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/HUVEC_LS_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed'))


##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 2000) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)


pdf('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/HUVEC_LS_VS_OS_specific_peaks_union_q0.05_fc_0.5_peaks_anno.pdf')
plotAnnoBar(peakAnnoList)
dev.off()


LS_speci <- as.data.frame(peakAnnoList$LS_speci)
OS_speci <- as.data.frame(peakAnnoList$OS_speci)

write.csv(LS_speci, file = "output.csv", quote = FALSE,  row.names = FALSE)








###########specific_peaks_overlaped_ATAC

files <- list(LS_speci_ATAC = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/overlaped_with_ATAC/LSS_DiffBind_DEseq2_specific_peaks_overlaped_ATAC_peaks.bed') , 
              OS_speci_ATAC = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/overlaped_with_ATAC/OSS_DiffBind_DEseq2_specific_peaks_overlaped_ATAC_peaks.bed'))


##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 2000) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)


pdf('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/overlaped_with_ATAC/HUVEC_LS_VS_OS_specific_peaks_overlaped_ATAC_union_q0.05_fc_0.5_peaks_anno.pdf')
plotAnnoBar(peakAnnoList)
dev.off()


LS_speci <- as.data.frame(peakAnnoList$LS_speci)
OS_speci <- as.data.frame(peakAnnoList$OS_speci)

write.csv(LS_speci, file = "H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/overlaped_with_ATAC/HUVEC_LS_specific_peaks_DiffBind_DESeq2_q0.05_fc0.5_overlaped_ATAC_Anno.csv", quote = TRUE,  row.names = FALSE)
write.csv(OS_speci, file = "H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/one-dimensional_new/DiffBind/Diffbind_DEseq2/overlaped_with_ATAC/HUVEC_OS_specific_peaks_DiffBind_DESeq2_q0.05_fc0.5_overlaped_ATAC_Anno.csv", quote = TRUE,  row.names = FALSE)










