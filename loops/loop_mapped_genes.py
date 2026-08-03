# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 16:03:04 2026

@author: lenovo
"""

import pandas as pd


s33_os_loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\IF_2+\\OS_s33_chr3_57950369_57953632_peaks_+-1kb_seq7858_rs268771_57952265_HUVEC_nonrisk_enhancer_speci.bedpe' , header = None)
s33_os_loops.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']



RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\HUVEC_DEGs_LS_VS_OS.csv' , header = 0)



chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']




RNA = RNA[RNA['Chr'].isin(chrom)]


DEGs = RNA[(RNA['padj'] <= 0.05) & ((RNA['log2FoldChange'] >= 0.5) | (RNA['log2FoldChange'] <= -0.5))]


LS_up = DEGs[DEGs['log2FoldChange'] >= 0.5]
LS_down = DEGs[DEGs['log2FoldChange'] <= -0.5]








