#!/bin/bash

outdir=/scratch/2023-12-18/bio-shenw/Ljniu/K562/K562_H/one-dimensional_2
rawdata_dir=/scratch/2023-12-18/bio-shenw/Ljniu/K562/K562_H/one-dimensional_2
genome_bowtie2_index=/scratch/2023-12-18/bio-shenw/ref/Human/hg38/bowtie2_index_bylxx/hg38.fa
threads=2

cd ${outdir}
for sample in ${rawdata_dir}/*_sorted.nodups.bam;do
    sample=$(basename $sample _sorted.nodups.bam)
    

    frac=$( samtools idxstats ${rawdata_dir}/${sample}_sorted.nodups.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=30000000/total; if (frac > 1) {print 1} else {print frac}}'  )
    echo "${frac}"
    bsub -J pairs -q ser -n ${threads} -o subset3000w.log -e subset3000w.err -R "span[ptile=20]" <<EOF
    #samtools view -O bam -s ${frac} ${rawdata_dir}/${sample}_sorted.nodups.bam > ${sample}.nodups_sorted_subset3000w.bam
    #samtools index ${sample}.nodups_sorted_subset3000w.bam
    samtools bam2fq  ${sample}.nodups_sorted_subset3000w.bam > ${sample}.nodups_sorted_subset3000w.fastq
    bowtie2  -p ${threads} --no-unal -x $genome_bowtie2_index -U ${sample}.nodups_sorted_subset3000w.fastq -S ${sample}_subset3000w_bowtie2.sam
    wait
    samtools sort -O bam -@ ${threads} -o ${sample}_sorted_subset3000w_bowtie2.bam ${sample}_subset3000w_bowtie2.sam
    wait
    samtools index ${sample}_sorted_subset3000w_bowtie2.bam
    wait

    
    #/work/bio-tanyj/miniconda3/bin/bamCoverage -b ${sample}.nodups_sorted_subset3000w.bam -o ${sample}_RPKM_10bp.bw -p 10 --binSize 10 --extendReads --minMappingQuality 20 --ignoreDuplicates --normalizeUsing RPKM
    macs2 callpeak -t ${sample}_sorted_subset3000w_bowtie2.bam -f BAM -g hs -n ${sample}_q0.05 -q 0.05 -B
    
    

EOF
done



