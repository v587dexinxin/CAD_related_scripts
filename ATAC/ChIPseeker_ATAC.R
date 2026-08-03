library(ChIPseeker)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("org.Mm.eg.db")
library("GenomicFeatures")

spompe <- makeTxDbFromGFF('H:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')
files <- list(WT = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/ATAC-seq/peaks/union_peaks/HUVEC_ATAC_WT_peaks_sorted_merged.bed') , 
              LS = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/ATAC-seq/peaks/union_peaks/HUVEC_ATAC_LS_peaks_sorted_merged.bed') , 
              OS = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/ATAC-seq/peaks/union_peaks/HUVEC_ATAC_OS_peaks_sorted_merged.bed'))


##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 2000) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)


pdf('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/HiRPC/ATAC-seq/peaks/union_peaks/HUVEC_WT_VS_LS_OS_ATAC_union_q0.05_peaks_anno.pdf')
plotAnnoBar(peakAnnoList)
dev.off()




