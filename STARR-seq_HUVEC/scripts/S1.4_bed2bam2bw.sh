#!/bin/bash
# variables need to be defined firstly, otherwise would run failed (if defined in the bsub workflow).

outdir=/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/workspace_nonrisk/mapping/bw
rawdata_dir=/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/workspace_nonrisk/mapping/UMI_dedup
referenceChrom=/scratch/2025-09-01/bio-shenw/ref/Human/hg38/hg38.chrom.size
NumThread=20

cd $outdir
for sample in ${rawdata_dir}/Uniq_UMI_NodT*.bed;do
	sample=$(basename $sample .bed)
	R1=$rawdata_dir/${sample}.bed
    R1_sorted=$rawdata_dir/${sample}_sorted.bed


	bsub -J bed2bam2bw -q ser -n 20 -R "span[ptile=20]" -o ${sample}.log -e ${sample}.err <<EOF
	echo -e "$(date '+%Y-%m-%d %H:%M:%S') $sample processing ... "
    
#1.sorted
    awk '{for(i=1;i<=NF;i++){if(\$i ~ /[eE]/){printf "%f\t", \$i}else{printf "%s\t", \$i}}; print ""}' ${R1} > $rawdata_dir/${sample}_1.bed
    sort -k1,1 -k2,2n -k3,3n  $rawdata_dir/${sample}_1.bed> ${R1_sorted}

#2. do bed2bw
    bedtools bedtobam -i ${R1_sorted} -g ${referenceChrom} > ${outdir}/${sample}.bam
    samtools  sort -@ ${NumThread} ${outdir}/${sample}.bam -o ${outdir}/${sample}_sorted.bam
    samtools index ${outdir}/${sample}_sorted.bam
    /work/bio-tanyj/miniconda3/bin/bamCoverage --normalizeUsing RPKM -p ${NumThread} -bs 10 -e 0 --minMappingQuality 20 -b ${outdir}/${sample}_sorted.bam -o ${outdir}/${sample}_RPKM_10bp.bw
    /work/bio-tanyj/miniconda3/bin/bamCoverage -p ${NumThread} -bs 10 -e 0 --minMappingQuality 20 -b ${outdir}/${sample}_sorted.bam -o ${outdir}/${sample}_10bp.bw

    macs2 callpeak -t ${outdir}/${sample}_sorted.bam  -f BAM -g hs -n ${sample}_q0.05 -q 0.05 -B


 	echo -e "$(date '+%Y-%m-%d %H:%M:%S') $sample bed2bam2bw ..."

	

EOF
done
 
