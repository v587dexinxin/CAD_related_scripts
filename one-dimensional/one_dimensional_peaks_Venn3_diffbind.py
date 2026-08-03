# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:13:40 2026

@author: lenovo
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3
import pandas as pd
import os
import pandas as pd 
import matplotlib
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3


def Load_bed(input_file):
    
    data = pd.read_table(input_file , header=None)
    data.columns = ['chr' , 'start' , 'end']
    return data


def Common_peaks(peaks1 , peaks2):
    n = 0 ; common = []
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    
    for g in chrom:
        tmp1 = peaks1[peaks1['chr'] == g]
        tmp2 = peaks2[peaks2['chr'] == g]
        for i in tmp1.index:
            start = tmp1.loc[i]['start']
            end = tmp1.loc[i]['end']
            mask = (tmp2['start'] <= end) & (tmp2['end'] >= start)
            overlap = tmp2[mask]
            if len(overlap) != 0:
                n += 1
                c = tuple(pd.concat([tmp1.loc[i] , overlap.iloc[0]] , axis = 0))
                common.append(c)
    print (n)
    common = pd.DataFrame(common)
    
    return common





def specific_peaks2(peaks1 , peaks2):
    n = 0 ; speci1 = pd.DataFrame([])
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    
    for g in chrom:
        tmp1 = peaks1[peaks1['chr'] == g]
        tmp2 = peaks2[peaks2['chr'] == g]
        for i in tmp1.index:
            start = tmp1.loc[i]['start']
            end = tmp1.loc[i]['end']
            mask = (tmp2['start'] <= end) & (tmp2['end'] >= start)
            overlap = tmp2[mask]
            if len(overlap) == 0:
                speci1 = pd.concat([speci1 , tmp1.loc[i:i]])
                n += 1
            
    print (n)
    
    return speci1


def plot_venn3(A_s , B_s , AB , C_s , AC , BC , ABC , title , out , lab1 , lab2 , lab3):
    fig = plt.figure(figsize = (10, 10))
    venn3(subsets=(A_s , B_s , AB , C_s , AC , BC , ABC), set_labels=(lab1, lab2 , lab3))
    for text_obj in fig.findobj(matplotlib.text.Text):
        text_obj.set_fontsize(20)
        
    plt.title(title , fontsize = 20)
    run_Plot(fig , out)
    
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    

WT = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_WT_allreps_q0.05_peaks_sorted_merged.bed')    
LS = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_LS_allreps_q0.05_peaks_sorted_merged.bed')    
OS = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_OS_allreps_q0.05_peaks_sorted_merged.bed')    


LSS_VS_OSS_LSS_s = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed')
LSS_VS_OSS_OSS_s = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed')
LSS_VS_OSS_common = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_common_peaks_q0.05_fc0.5_clean.bed')



WT_VS_OSS_WT_s = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\WT_VS_OS\\HUVEC_WT_VS_OS_DESeq2_WT_specific_peaks_q0.05_fc0.5_clean.bed')
WT_VS_OSS_OS_s = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\WT_VS_OS\\HUVEC_WT_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed')
WT_VS_OSS_common = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\WT_VS_OS\\HUVEC_WT_VS_OS_DESeq2_common_peaks_q0.05_fc0.5_clean.bed')




WT_VS_LSS_WT_s = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\WT_VS_LS\\HUVEC_WT_VS_LS_DESeq2_WT_specific_peaks_q0.05_fc0.5_clean.bed')
WT_VS_LSS_LS_s = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\WT_VS_LS\\HUVEC_WT_VS_LS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed')
WT_VS_LSS_common = Load_bed('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\WT_VS_LS\\HUVEC_WT_VS_LS_DESeq2_common_peaks_q0.05_fc0.5_clean.bed')


ABC = 34708

ABc = len(Common_peaks(WT_VS_OSS_WT_s , LSS_VS_OSS_LSS_s))
AbC = len(Common_peaks(WT_VS_LSS_WT_s , LSS_VS_OSS_OSS_s))
aBC = len(Common_peaks(WT_VS_LSS_LS_s , WT_VS_OSS_OS_s))



Abc = len(WT) - ABC - ABc - AbC
aBc = len(LS) - ABC - ABc - aBC
abC = len(OS) - ABC - aBC - AbC




out = 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\one-dimensional\\WT_VS_LS_VS_OS_Hi-Coatis_peaks_DiffBind_Venn3.pdf'
plot_venn3(Abc , aBc , ABc , abC , AbC , aBC , ABC , 'WT_VS_LS_VS_OS' , out , 'WT' , 'LS' , 'OS')







