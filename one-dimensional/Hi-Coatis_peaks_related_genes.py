# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 15:37:44 2026

@author: lenovo
"""

import pandas as pd


def peak_related_genes(peaks):
    peak_genes = []
    for g in chrom:
        tmp_peaks = peaks[peaks['chr'] == g]
        tmp_rna = RNA[RNA['Chr'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            mask = (start <= tmp_rna['pro_e']) & (end >= tmp_rna['pro_s'])
            overlap = tmp_rna[mask]
            if len(overlap) > 0:
                genes = overlap.iloc[0]
                g_s = genes['Start']
                g_e = genes['End']
                g_name = genes['Gene_Name']
                g_score = '.'
                g_strand = genes['Strand']
                peak_genes.append((g , g_s , g_e , g_name , g_score , g_strand))
    peak_genes = pd.DataFrame(peak_genes , columns = ['chr' , 'start' , 'end' , 'Gene_Name' , 'score' , 'strand'])
    peak_genes = peak_genes.drop_duplicates(subset=["Gene_Name"])
    
    return peak_genes






lss_speci = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
lss_speci.columns = ['chr' , 'start' , 'end']
oss_speci = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
oss_speci.columns = ['chr' , 'start' , 'end']


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




chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


lss_peak_genes = peak_related_genes(lss_speci)
oss_peak_genes = peak_related_genes(oss_speci)


lss_peak_genes.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\specific_peaks_related_genes\\LSS_specific_peaks_related_gene.bed' , header = None , index = None , sep = '\t')
oss_peak_genes.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\specific_peaks_related_genes\\OSS_specific_peaks_related_gene.bed' , header = None , index = None , sep = '\t')




















