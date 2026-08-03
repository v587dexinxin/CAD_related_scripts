# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 16:07:23 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()



CAD_SNPs = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_SNPs_LD0.99_all.bed' , header = 0)

peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_WT_allreps_q0.05_peaks_sorted_merged.narrowPeak' , header = None , usecols = (0 , 1 , 2))
peaks.columns = ['chr' , 'start' , 'end']
peaks_562 = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak' , header = None , usecols = (0 , 1 , 2))
peaks_562.columns = ['chr' , 'start' , 'end']
peaks_116 = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_WT_allreps_q0.05_peaks_sorted_merged.narrowPeak' , header = None , usecols = (0 , 1 , 2))
peaks_116.columns = ['chr' , 'start' , 'end']
peaks_H3K27ac = pd.read_table('H:\\work\\literature_data\\HUVEC\\HUVEC_HeK27ac_hg38_ENCFF077LGZ.bed' , header = None , usecols = (0 , 1 , 2))
peaks_H3K27ac.columns = ['chr' , 'start' , 'end']
peaks_H3K4me3 = pd.read_table('H:\\work\\literature_data\\HUVEC\\HUVEC_H3K4me3_hg38_ENCFF550OWZ.bed' , header = None , usecols = (0 , 1 , 2))
peaks_H3K4me3.columns = ['chr' , 'start' , 'end']
peaks_ATAC = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\HUVEC_ATAC_union_peaks_sorted_merged.narrowPeak' , header = None , usecols = (0 , 1 , 2)) 
peaks_ATAC.columns = ['chr' , 'start' , 'end']


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

n = 0
for g in chrom:
    print (g)
    tmp_snps = CAD_SNPs[CAD_SNPs['CHR_ID'] == g.lstrip('chr')]
    tmp_peaks = peaks[peaks['chr'] == g]
    tmp_peaks_562 = peaks_562[peaks_562['chr'] == g]
    tmp_peaks_116 = peaks_116[peaks_116['chr'] == g]
    for i in tmp_snps.index:
        pos = tmp_snps.loc[i]['CHR_POS']
        pos_s = pos - 1000
        pos_e = pos + 1000
        mask = (tmp_peaks['start'] <= pos_e) & (tmp_peaks['end'] >= pos_s)
        mask1 = (tmp_peaks_562['start'] <= pos_e) & (tmp_peaks_562['end'] >= pos_s)
        mask2 = (tmp_peaks_116['start'] <= pos_e) & (tmp_peaks_116['end'] >= pos_s)
        overlap = tmp_peaks[mask]
        overlap1 = tmp_peaks_562[mask1]
        overlap2 = tmp_peaks_116[mask2]
        if (len(overlap) > 0) or (len(overlap1) > 0) or (len(overlap2) > 0):
            n += 1







m = 0
for g in chrom:
    print (g)
    tmp_snps = CAD_SNPs[CAD_SNPs['CHR_ID'] == g.lstrip('chr')]
    tmp_peaks = peaks[peaks['chr'] == g]
    tmp_peaks_562 = peaks_562[peaks_562['chr'] == g]
    tmp_peaks_116 = peaks_116[peaks_116['chr'] == g]
    tmp_peaks_H3K27ac =  peaks_H3K27ac[peaks_H3K27ac['chr'] == g]
    tmp_peaks_H3K4me3 =  peaks_H3K4me3[peaks_H3K4me3['chr'] == g]
    tmp_peaks_atac = peaks_ATAC[peaks_ATAC['chr'] == g]
    
    for i in tmp_snps.index:
        pos = tmp_snps.loc[i]['CHR_POS']
        pos_s = pos - 1000
        pos_e = pos + 1000
        mask = (tmp_peaks['start'] <= pos_e) & (tmp_peaks['end'] >= pos_s)
        mask1 = (tmp_peaks_562['start'] <= pos_e) & (tmp_peaks_562['end'] >= pos_s)
        mask2 = (tmp_peaks_116['start'] <= pos_e) & (tmp_peaks_116['end'] >= pos_s)
        mask3 = (tmp_peaks_H3K27ac['start'] <= pos_e) & (tmp_peaks_H3K27ac['end'] >= pos_s)
        mask4 = (tmp_peaks_H3K4me3['start'] <= pos_e) & (tmp_peaks_H3K4me3['end'] >= pos_s)
        mask5 = (tmp_peaks_atac['start'] <= pos_e) & (tmp_peaks_atac['end'] >= pos_s)
        
        overlap = tmp_peaks[mask]
        overlap1 = tmp_peaks_562[mask1]
        overlap2 = tmp_peaks_116[mask2]
        overlap3 = tmp_peaks_H3K27ac[mask3]
        overlap4 = tmp_peaks_H3K4me3[mask4]
        overlap5 = tmp_peaks_atac[mask5]
        
        if (len(overlap) > 0) or (len(overlap1) > 0) or (len(overlap2) > 0) or (len(overlap3) > 0) or (len(overlap4) > 0) or (len(overlap5) > 0):
            m += 1






# 数据
labels = ['CAD SNPs']
in_peaks = [m]
not_in_peaks = [len(CAD_SNPs) - m]

# 颜色
color_peak = '#FFD700'      # 亮黄色
color_not_peak = '#1F4E79'  # 深蓝色

fig, ax = plt.subplots(figsize=(4,6))

# 画堆叠柱状图
bar1 = ax.bar(labels, in_peaks, color=color_peak, label='SNP in Hi-Coatis peaks')
bar2 = ax.bar(labels, not_in_peaks, bottom=in_peaks, color=color_not_peak,
              label='SNP not in Hi-Coatis peaks')

# 添加数字
ax.text(0, in_peaks[0]/2, str(in_peaks[0]), ha='center', va='center', fontsize=12)
ax.text(0, in_peaks[0] + not_in_peaks[0]/2, str(not_in_peaks[0]),
        ha='center', va='center', color='white', fontsize=12)

# 坐标轴
ax.set_ylabel('Number of SNPs')
ax.set_title('Distribution of CAD-associated SNPs')

# 图例
ax.legend()

plt.tight_layout()
plt.show()







