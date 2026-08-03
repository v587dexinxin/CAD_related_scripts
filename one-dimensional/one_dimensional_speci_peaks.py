# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:13:42 2026

@author: lenovo
"""

import pandas as pd 


deseq2_pos = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header=0)

# diff = deseq2_pos[deseq2_pos['padj'] <= 0.05]

# oss_speci = diff[diff['log2FoldChange'] >= 0.5]
# lss_speci = diff[diff['log2FoldChange'] < -0.5]



no_sign = deseq2_pos[deseq2_pos['direction'] == 'not_significant']
lss_speci = deseq2_pos[deseq2_pos['direction'] == 'LSS_specific']
oss_speci = deseq2_pos[deseq2_pos['direction'] == 'OSS_specific']




no_sign.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_common_peaks_q0.05_fc0.5.csv' , index = None)
no_sign[['chr' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_common_peaks_q0.05_fc0.5.bed' , header = None , index = None , sep = '\t')













