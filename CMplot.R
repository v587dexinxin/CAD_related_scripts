library("CMplot")
data = read.table('D:/work/Postdoctoral/GWAS¼²²¡Î»µã¼ì²â/results/CAD/Hapmap/CAD_related_SNPs_LD0.9_all.bed')
CMplot(data,plot.type="d",bin.size=2000,chr.den.col=c("darkgreen", "red", "yellow"),file="pdf",memo="GM_SNP_2K", file.output=TRUE,verbose=TRUE, bin.range = c(1 , 3))





