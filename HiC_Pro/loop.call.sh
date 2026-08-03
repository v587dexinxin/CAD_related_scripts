#!/bin/bash
PROGRAM_DIRECTORY='/work/bio-shenw/tools/ChIAPET-callloops_huangxy'
OUTPUT_DIRECTORY='/scratch/2024-02-26/bio-shenw/Ljniu/K562/K562_H_New/loops/ChIA-pipline_8K+'
OUTPUT_PREFIX='K562_Hemin_New_merged4.hg38.nodups.pairs_+8K'

EXTENSION_LENGTH='500'
CHROM_SIZE_INFO='/scratch/2024-02-26/bio-shenw/ref/Human/hg38/hg38.chrom.size'


#java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.file.Pet2Cluster1 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt ${EXTENSION_LENGTH} ${CHROM_SIZE_INFO}

java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.file.Pet2Cluster1 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt ${EXTENSION_LENGTH} ${CHROM_SIZE_INFO}
LANG=C sort -k1,1 -k4,4  < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.txt >  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted
java -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.chiapet.PetCluster2 ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pre_cluster.sorted ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster

awk '{if($13>=2)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered
cut -f1-3,7-9,13-15 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.txt

###### calculation of p-values
cut -f1-3 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln
cut -f4-6 < ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.ipet >>${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln

cut -f1-3,13 <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1
cut -f7-9,13 <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2
for y in anchor1 anchor2
do
    java -Xmx10G -cp ${PROGRAM_DIRECTORY}/LGL.jar LGL.shortReads.TagCountInGivenRegions ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.${y} ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.${y}.tagCount.txt 1 2
done

wc -l ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.aln | sed 's/ /\t/g' |  cut -f1 >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt
## calculate p-value
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.nTags.txt nTags.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor1.tagCount.txt  ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered.anchor2.tagCount.txt |  cut -f4,5,10  >${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt
cp ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt data.txt
R --vanilla < ${PROGRAM_DIRECTORY}/hypergeometric.r
mv -f result.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt
paste ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.filtered ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.petCount.tagCount.txt ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.pvalue.hypergeo.txt | cut -f1-3,7-9,13-15,19-24 > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt
awk -v pvalue_cutoff=${PVALUE_CUTOFF_INTERACTION=0.05} '{if($13<pvalue_cutoff+0.0)print}' <${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.withpvalue.txt > ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt
date '+%s' >> ${OUTPUT_DIRECTORY}/time.txt
awk '{print $7"\t"$8}' ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.cluster.FDRfiltered.txt |LANG=C sort |uniq -c | awk '{if($2>=10){a+=$1;b+=$1*$3;f+=$1;g+=$1*$3}else{c[$2]+=$1;d[$2]+=$1*$3;f+=$1;g+=$1*$3}}END{for(i in c){print c[i]"\t"i"\t"d[i]};print a"\t10\t"b;print f"\t11\t"g}' |LANG=C sort -k2,2n | awk 'BEGIN{print "PET counts\tNo. of clusters\tNo.intra-chrom clusters\tNo.inter-chrom clusters\tPercent of intra-chrom clusters"}{if($2==10){print ">="$2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else if($2==11){print "Total\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}else{print $2"\t"$1"\t"$3"\t"$1-$3"\t"$3/$1}}'> ${OUTPUT_DIRECTORY}/${OUTPUT_PREFIX}.PET_count_distribution.txt



