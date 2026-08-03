# -*- coding: utf-8 -*-
"""
Created on Thu May 14 19:32:01 2026

@author: lenovo
"""

import pandas as pd



housekeeping = ['GAPDH' , 'ACTB' , 'B2M' , 'HPRT1' , 'TBP' , 'RPLP0' , 'RPL13A' , 'YWHAZ' , 'PPIA' , 'UBC']

RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)
union_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\LS_OS_union_peaks_sorted_merged.bed' , header= None)
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
        
        








h_pro.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\LS_OS_union_peaks_sorted_merged_housekeeping.bed' , header = None , index = None , sep = '\t')



















