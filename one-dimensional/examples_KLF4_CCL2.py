# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:05:36 2026

@author: lenovo
"""



import pandas as pd


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


##----------------RNA---------------------

RNA = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0 , sep = ',')
RNA = RNA.drop_duplicates(subset = ['Gene_Name'] , keep = 'first')
RNA = RNA[~RNA["Gene_Name"].str.contains("STRG", na=False)]
RNA['f_s'] = RNA['Start'] - 2000
RNA['f_e'] = RNA['End'] + 2000

a = RNA[RNA['Strand'] == '+']
a['pro_s'] = a['Start'] - 2000
a['pro_e'] = a['Start'] + 500

b = RNA[RNA['Strand'] == '-']
b['pro_s'] = b['End'] - 500
b['pro_e'] = b['End'] + 2000

RNA = pd.concat([a , b])

RNA['LS_FPKM'] = (RNA['LS_R1_FPKM'] + RNA['LS_R2_FPKM']) / 2
RNA['OS_FPKM'] = (RNA['OS_R1_FPKM'] + RNA['OS_R2_FPKM']) / 2

expressed_rna = RNA[(RNA['LS_FPKM'] >= 2) | (RNA['OS_FPKM'] >= 2)]





klf4 = RNA[RNA['Gene_Name'] == 'KLF4']
ccl2 = RNA[RNA['Gene_Name'] == 'CCL2']




##----------------Hi-Coatis---------------------


count = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DiffBind_raw_counts_for_DESeq2.csv' , header = 0)
count = count[count['CHR'].isin(chrom)]


selected_peaks = [('chr9' , 107627730 , 107628066 , 'KLF4') , ('chr17' , 34382438 , 34383171 , 'CCL2')]

selected_count = []

for i in selected_peaks:
    g = i[0]
    start = i[1]
    end = i[2]
    gene_name = i[3]
    tmp_count = count[count['CHR'] == g]
    mask = (tmp_count['START'] <= end) & (tmp_count['END'] >= start)
    overlap = tmp_count[mask]
    if len(overlap) == 1:
        selected_count.append((overlap.iloc[0]['LS_R1'] , overlap.iloc[0]['LS_R2'] , overlap.iloc[0]['LS_R3'] , overlap.iloc[0]['OS_R1'] , overlap.iloc[0]['OS_R2'] , overlap.iloc[0]['OS_R3']))
        
        
        
selected_count = pd.DataFrame(selected_count)
        
        
        
        
        
##----------------ATAC-seq---------------------


count = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\union_peaks_related_ATAC\\HUVEC_ATAC_LS_VS_OS_DiffBind_raw_counts_for_DESeq2_Hi-Coatis_peaks.csv' , header = 0)
count = count[count['CHR'].isin(chrom)]


selected_peaks = [('chr9' , 107627730 , 107628066 , 'KLF4') , ('chr17' , 34382438 , 34383171 , 'CCL2')]

selected_count = []

for i in selected_peaks:
    g = i[0]
    start = i[1]
    end = i[2]
    gene_name = i[3]
    tmp_count = count[count['CHR'] == g]
    mask = (tmp_count['START'] <= end) & (tmp_count['END'] >= start)
    overlap = tmp_count[mask]
    if len(overlap) == 1:
        selected_count.append((overlap.iloc[0]['LS_R1'] , overlap.iloc[0]['LS_R2'] , overlap.iloc[0]['OS_R1'] , overlap.iloc[0]['OS_R2']))
        
        
        
selected_count = pd.DataFrame(selected_count)
        
        
        
        
        






































































