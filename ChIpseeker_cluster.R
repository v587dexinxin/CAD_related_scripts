library(ChIPseeker)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("org.Mm.eg.db")
library("GenomicFeatures")

spompe <- makeTxDbFromGFF('D:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')
files = list(Cardiovascular = c('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99+-1bp_all.bed'))



##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 0) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)

pdf('D:/work/Postdoctoral/union_peakAnno_first6000.pdf')
plotAnnoBar(peakAnnoList)
dev.off()



peakAnnoList[["Cardiovascular"]]



cluster1 <- as.data.frame(peakAnnoList[["Cardiovascular"]]@anno)
c1 <- with(cluster1,grep("Promoter",annotation))
c1_intergenic <- cluster1[c1,]

write.table(c1_intergenic , "D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Promoter_1.txt" , sep="\t" , quote=F , row.names=F)


write.table(c1_intergenic[,c('seqnames' , 'start' , 'end')] , "D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Promoter_1.bed" , sep="\t" , quote=F , row.names=F , col.names = F)



c2 <- with(cluster1,grep("Exon",annotation))
c2_intergenic <- cluster1[c2,]

write.table(c2_intergenic , "D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Exon_1.txt" , sep="\t" , quote=F , row.names=F)


write.table(c2_intergenic[,c('seqnames' , 'start' , 'end')] , "D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Exon_1.bed" , sep="\t" , quote=F , row.names=F , col.names = F)


c3 <- with(cluster1,grep("Intron",annotation))
c3_intergenic <- cluster1[c3,]
c4 <- with(cluster1,grep("Downstream",annotation))
c4_intergenic <- cluster1[c4,]
c5 <- with(cluster1,grep("Distal Intergenic",annotation))
c5_intergenic <- cluster1[c5,]

c345_intergenic<- rbind(c3_intergenic , c4_intergenic , c5_intergenic)


write.table(c345_intergenic , "D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Noncoding_1.txt" , sep="\t" , quote=F , row.names=F)


write.table(c345_intergenic[,c('seqnames' , 'start' , 'end')] , "D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Noncoding_1.bed" , sep="\t" , quote=F , row.names=F , col.names = F)


