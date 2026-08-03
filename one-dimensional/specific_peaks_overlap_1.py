# -*- coding: utf-8 -*-
"""
Created on Fri May 22 12:49:09 2026

@author: lenovo
"""

import pandas as pd



data = pd.read_table('H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/one-dimensional_new/bam/DESeq2_only/HiCoatis_LSS_vs_OSS_DESeq2_results.txt' , header=0)



diff = data[data['padj'] <= 0.05]

lss_speci = diff[diff['log2FoldChange'] >= 1]
oss_speci = diff[diff['log2FoldChange'] < -1]



oss_speci[['Chr' , 'Start' , 'End']].to_csv('H:/work/Postdoctoral/GWAS疾病位点检测/results/HiRPC/one-dimensional_new/bam/DESeq2_only/HiCoatis_LSS_vs_OSS_DESeq2_OSS_speci.bed' , header = None , index = None , sep = '\t')






deseq2_pos = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header=0)

deseq2_pos.columns = ['peak_id', 'Chr', 'Start', 'End', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'direction']

diff = deseq2_pos[deseq2_pos['padj'] <= 0.05]

oss_speci_1 = diff[diff['log2FoldChange'] >= 0.5]
lss_speci_1 = diff[diff['log2FoldChange'] < -0.5]




speci1 = lss_speci
speci2 = lss_speci_1


speci1 = oss_speci_1
speci2 = oss_speci



oss_edgR_speci = []
common = []

n = 0
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
for g in chrom:
    tmp1 = speci1[speci1['Chr'] == g]
    tmp2 = speci2[speci2['Chr'] == g]
    for i in tmp1.index:
        start = tmp1.loc[i]['Start']
        end = tmp1.loc[i]['End']
        mask = (tmp2['Start'] <= end) & (tmp2['End'] >= start)
        overlap = tmp2[mask]
        if len(overlap) > 0:
            n += 1
            common.append((g , start , end))
        else:
            oss_edgR_speci.append((g , start , end))
            
oss_edgR_speci = pd.DataFrame(oss_edgR_speci)
common = pd.DataFrame(common)


print (len(speci1) - n , n , len(speci2) - n)






lss_speci[['Chr', 'Start', 'End']].to_csv()











