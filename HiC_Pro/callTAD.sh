#!/bin/bash
#BSUB -J TADCall
#BSUB -e TADCall.e
#BSUB -o TADCall.o
#BSUB -q ser 
#BSUB -n 12
#BSUB -R "span[ptile=12]"
source /work/med-wulf/.bashrc
cd /scratch/2024-05-06/med-wulf/Hi-C-Result/
path="/scratch/2024-05-06/med-wulf/Hi-C-Result/raw_data"
for name in "KHM-2" "KHM-WT";do
/data/med-wulf/software/HiC-Pro/build/HiC-Pro-master/bin/utils/hicpro2juicebox.sh -i /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/hic_results/data/${name}/${name}.allValidPairs -g /data/med-wulf/Reference/hg38.chrom.sizes -j /data/med-wulf/software/juicer/juicer_tools_1.22.01.jar -r /data/med-wulf/Reference/hg38_DpnII.bed -o /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/
#add normalization data in hic file.
java -jar /data/med-wulf/software/juicer/juicer_tools_1.22.01.jar addNorm -w 10000 -k KR /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs.hic
done
for name in "KHM-2" "KHM-WT";do
	#java -jar /data/med-wulf/software/juicer/juicer_tools_1.22.01.jar eigenvector KR /data/med-wulf/WZ/Hi-C-Result/HiC/${name}.allValidPairs.hic chr1 BP 1000000 /data/med-wulf/WZ/Hi-C-Result/HiC/${name}_1000000.eigen.txt
	#TAD, for the contact map is too sparse, so need to use --ignore-sparsity (but we need to merge replicates to get high resolution)
	for res in "10000" "25000";do
	java -jar /data/med-wulf/software/juicer/juicer_tools_1.22.01.jar arrowhead -k KR --threads 12 -r ${res} --ignore-sparsity /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs.hic /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}_${res}_arrowhead
	done
	#loop, for the contact map is too sparse, so need to use --ignore-sparsity (but we need to merge replicates to get high resolution)
	java -jar /data/med-wulf/software/juicer/juicer_tools_1.22.01.jar hiccups --cpu -m 1024 -r 5000,10000 --ignore-sparsity --threads 12 /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs.hic /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}_hiccups_loops
done
#due to low depth, change to hicexplorer to call TAD and loops.
conda activate hicexplorer
for name in "KHM-2" "KHM-WT";do
hicConvertFormat -m /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs.hic --inputFormat hic --outputFormat cool -o /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs.cool --resolutions 10000
hicCorrectMatrix correct --matrix /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs_10000.cool --correctionMethod KR --outFileName /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}_corrected.cool
hicFindTADs -m /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}_corrected.cool --correctForMultipleTesting fdr --outPrefix /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name} --numberOfProcessors 8
hicDetectLoops -m /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.allValidPairs_10000.cool -o /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05
done
#calculate the number of loops and TADs in each sample.
rm /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/24.5.2.hic.statistic.txt
for name in "KHM-2" "KHM-WT";do
	TADNum=$(cat /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}_domains.bed | wc -l)
	LoopNum=$(cat /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/${name}.bedgraph | wc -l)
	echo $name $TADNum $LoopNum >> /scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/HiC/24.5.2.hic.statistic.txt
done

