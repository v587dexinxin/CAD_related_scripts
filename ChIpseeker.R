library(ChIPseeker)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("org.Mm.eg.db")
library("GenomicFeatures")

spompe <- makeTxDbFromGFF('D:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')
files = list(Cardiovascular = c('D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/Hapmap/CAD_related_SNPs_LD0.9+-1bp_all.bed'))



##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 0) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,genomicAnnotationPriority = c("Promoter", "Exon", "Intron","Intergenic"),verbose = FALSE)

pdf('D:/work/Postdoctoral/union_peakAnno_5.pdf')
plotAnnoBar(peakAnnoList)
dev.off()



peakAnnoList[["Cardiovascular"]]



cluster1 <- as.data.frame(peakAnnoList[["Cardiovascular"]]@anno)
c1 <- with(cluster1,grep("Promoter",annotation))
c1_intergenic <- cluster1[c1,]

write.table(c1_intergenic , "D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_related_SNPs_Promoter.txt" , sep="\t" , quote=F , row.names=F)


write.table(c1_intergenic[,c('seqnames' , 'start' , 'end')] , "D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_related_SNPs_Promoter.bed" , sep="\t" , quote=F , row.names=F , col.names = F)

pro_file = list(Promoter = c('D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_related_SNPs_Promoter_+-60bp.bed'))

plotAvgProf2(pro_file, TxDb=spompe, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)
			 
			 
			 
####promoter			 
			 
spompe <- makeTxDbFromGFF('D:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')
files = list(Cardiovascular = c('D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_disease/promoter/Cardiovascular_promoter_related_SNPs_LD0.8+-1bp_all.bed'))



##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 0) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,genomicAnnotationPriority = c("Promoter", "Exon", "Intron","Intergenic"),verbose = FALSE)

pdf('D:/work/Postdoctoral/union_peakAnno_5.pdf')
plotAnnoBar(peakAnnoList)
dev.off()



peakAnnoList[["Cardiovascular"]]



cluster1 <- as.data.frame(peakAnnoList[["Cardiovascular"]]@anno)
c1 <- with(cluster1,grep("Promoter",annotation))
c1_intergenic <- cluster1[c1,]

write.table(c1_intergenic , "D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_disease/promoter/Cardiovascular_promoter_related_SNPs_Promoter_LD0.8.txt" , sep="\t" , quote=F , row.names=F)


write.table(c1_intergenic[,c('seqnames' , 'start' , 'end')] , "D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_disease/promoter/Cardiovascular_promoter_related_SNPs_Promoter_LD0.8.bed" , sep="\t" , quote=F , row.names=F , col.names = F)

tmp <- c1_intergenic[,c('seqnames' , 'start' , 'end')]
x <- c(60 , -60)
tmp[,c('start' , 'end')] <- sweep(tmp[,c('start' , 'end')], 2, x)
write.table(tmp , "D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_disease/promoter/Cardiovascular_promoter_related_SNPs_Promoter_LD0.8_+-60bp.bed" , sep="\t" , quote=F , row.names=F , col.names = F)




pro_file = list(Promoter = c('D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/Cardiovascular_disease/promoter/Cardiovascular_promoter_related_SNPs_Promoter_LD0.8_+-60bp.bed'))


pdf('D:/work/Postdoctoral/promoter_SNP_TSS_distance_LD0.8.pdf')

plotAvgProf2(pro_file, TxDb=spompe, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)
			 
			 
dev.off()