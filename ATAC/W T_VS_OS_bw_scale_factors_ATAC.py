# -*- coding: utf-8 -*-
"""
Created on Mon May 18 18:50:05 2026

@author: lenovo
"""


import pandas as pd
import numpy as np




#######stable_peaks

# deseq2_pos = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\DiffBind\\DiffBind_DEseq2\\HUVEC_LS_VS_OS_ATAC_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header=0)

# deseq2_pos.columns = ['peak_id', 'chr', 'start', 'end', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'direction']

# no_signific = deseq2_pos[deseq2_pos['direction'] == 'not_significant']


no_diff = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\DiffBind\\DiffBind_DEseq2\\WT_VS_OS\\HUVEC_WT_VS_OS_ATAC_DESeq2_common_peaks_p0.05_fc0.5.csv' , header=0)


########housekeeping

housekeeping = ['GAPDH' , 'ACTB' , 'B2M' , 'HPRT1' , 'TBP' , 'RPLP0' , 'RPL13A' , 'YWHAZ' , 'PPIA' , 'UBC']

RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)
union_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_LS_OS_union_peaks_sorted_merged.bed' , header= None)
union_peaks.columns = ['chr' , 'start' , 'end']



h_pro = pd.DataFrame([])
for gene_name in housekeeping:
    gene = RNA[RNA['Gene_Name'] == gene_name]
    g = gene.iloc[0]['Chr']
    start = gene.iloc[0]['Start']
    end = gene.iloc[0]['End']
    strand = gene.iloc[0]['Strand']
    if strand == '+':
        pro_s = start - 2000
        pro_e = start + 500
    else:
        pro_s = end - 500
        pro_e = end + 2000
    tmp_peaks = union_peaks[union_peaks['chr'] == g]
    mask = (tmp_peaks['start'] <= pro_e) & (tmp_peaks['end'] >= pro_s)
    overlap = tmp_peaks[mask]
    if len(overlap) > 0:
        h_pro = pd.concat([h_pro , overlap])
        
        
# chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

# for g in chrom:
#     tmp_h = h_pro[h_pro['chr'] == g]
#     tmp_no_diff = no_diff[no_diff['chr'] == g]
#     if len(tmp_h) > 0:
#         for i in tmp_h.index:
#             start = tmp_h.loc[i]['start']
#             end = tmp_h.loc[i]['end']
#             mask = (tmp_no_diff['start'] <= end) & (tmp_no_diff['end'] >= start)
#             overlap = tmp_no_diff[mask]
#             if len(overlap) > 0:
#                 print (g , start , end)

no_diff_h = pd.concat([no_diff[['chr' , 'start' , 'end']] , h_pro])



no_diff_h.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\signals\\WT_OSS_ATAC_stable_housekeeping_peaks.bed' , header = None , index = None , sep = '\t')


############################Linux order.txt
#sort -k1,1V -k2,2n  LSS_OSS_stable_housekeeping_peaks.bed >LSS_OSS_stable_housekeeping_peaks_sorted.bed
#bedtools merge -i LSS_OSS_stable_housekeeping_peaks_sorted.bed > LSS_OSS_stable_housekeeping_peaks_sorted_merged.bed
#awk 'BEGIN{OFS="\t"}$1 ~ /^chr[0-9XYM]+$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 >= 0 && $3 > $2 {print $1, int($2), int($3)}' LSS_OSS_ATAC_stable_housekeeping_peaks_sorted_merged.bed | sort -k1,1 -k2,2n > LSS_OSS_ATAC_stable_housekeeping_peaks_sorted_merged_clean.bed
###清洗非法行
#bsub -J stable -q ser -n 20 -e stable.err -o stable.log -R "span[hosts=1]" multiBigwigSummary BED-file --BED LSS_OSS_stable_housekeeping_peaks_sorted_merged.bed -b ../HiRPC_LS_allreps_RPKM_10bp.bw ../HiRPC_OS_allreps_RPKM_10bp.bw --outRawCounts stable_peaks_bw_signal.tab -out stable_peaks_bw_signal.npz -p 20




#############################stable peaks-based scale factor



import pandas as pd
import numpy as np
from scipy.stats import trim_mean

signal_file = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\signals\\stable_housekeeping_rescaled_bw\\WT_OS_ATAC_stable_peaks_bw_signal.tab"

df = pd.read_csv(signal_file, sep="\t", header = 0)

signal_cols = df.columns[3:]

df = df[df[signal_cols].sum(axis=1) > 0].copy()

trimmed_means = {}

for col in signal_cols:
    values = df[col].dropna().values
    values = values[values > 0]
    trimmed_means[col] = trim_mean(values, proportiontocut=0.1)

trimmed_means = pd.Series(trimmed_means)

target = trimmed_means.mean()
scale_factors = target / trimmed_means

print("Trimmed mean signal on stable peaks:")
print(trimmed_means)

print("\nStable peak-based scale factors:")
print(scale_factors)

scale_factors.to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\signals\\stable_housekeeping_rescaled_bw\\WT_OSS_stable_peaks_trimmedMean_scale_factors.txt",
    sep="\t",
    header=["scale_factor"]
)


####results: LSS_scales factor: 0.882

