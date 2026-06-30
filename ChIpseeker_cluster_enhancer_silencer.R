library(ChIPseeker)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("org.Mm.eg.db")
library("GenomicFeatures")

spompe <- makeTxDbFromGFF('H:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')

files <- list(all_SNPs = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/CAD_related_Vascular endothelium diseasesSNPs_LD0.99+-1bp_all.bed') ,
              enhancer_116 = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/DESeq2/enhancer_silencers/position/HCT116_enhancer_deseq2_pos_+-1bp.bed') , 
              enhancer_huvec = c('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/DESeq2/enhancer_silencers/position/HUVEC_enhancer_deseq2_pos_+-1bp.bed') )


##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 500) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)


pdf('H:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/first_6000/plots/union_peakAnno_first6000_enhancer_silencers.pdf')
plotAnnoBar(peakAnnoList)
dev.off()

