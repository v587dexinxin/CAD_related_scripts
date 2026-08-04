#!/bin/bash
#BSUB -J UMItab                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 1                 ### ask for number of cores (default: 1)
#BSUB -W 35:3             ### set walltime limit: hh:mm
#BSUB -e umi.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o umi.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[hosts=1]"    ### ask for 40 cores per node

cd /scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/data/nonrisk/UMI_fq
outdir="/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/workspace_nonrisk/mapping/UMItab"

# UMI table (Read2 8bp)
 # output a tsv file which contains two columns: readID and UMI sequence.

for i in `ls | grep "_R2"`
do
    {
    id=${i%_R2_nonrisk.fastq.gz}
    #echo "$id,$i"
    zcat ${i} | \
        awk 'NR%4==1 {readID=$1; sub("^[^:]+:[^:]+:[^:]+:[^:]+:", "", readID)} NR%4==2 {UMI=$0; print readID, UMI}' > ${outdir}/${id}.UMItab
    }&
done
wait

