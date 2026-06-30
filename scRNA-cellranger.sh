for idx in {mh018,mh020,mh044,mh046};
do
    cd /public/home/zhangzb/monkey-heart/log
    qsub << EOF
        #PBS -N ${idx}_cellranger
        #PBS -l nodes=comput2:ppn=5
        #PBS -l walltime=240:00:00
        cd /public/home/zhangzb/Referencegenome/monkey/ensembl/cellranger_fascicularis5.0/cellranger

        cellranger count --id=${idx} \
                         --transcriptome=/public/home/zhangzb/Referencegenome/monkey/ensembl/cellranger_fascicularis5.0/Macaca_fascicularis_genome \
                         --fastqs=/public/home/zhangzb/monkey-heart/sn_seq/seq-data/Cleandata/${idx} \
                         --sample=${idx}-1,${idx}-2 \
                         --localcores=5
EOF
sleep 1
done
         