# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 15:00:15 2026

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



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
def Load_peaks(file , peaks_type):
    sz = os.path.getsize(file)
    if sz != 0:
        peaks = pd.read_table(file , header = None)
        if peaks_type == 'narrow':
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue' , 'lengtn']
        else:
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue']
    else:
        peaks = []

    return peaks

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



def specific_peaks3(peaks1 , peaks2 , peaks3):
    n = 0 ; speci1 = pd.DataFrame([]) ; speci2 = pd.DataFrame([])
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
                
    for g in chrom:
        tmp3 = speci1[speci1['chr'] == g]
        tmp4 = peaks3[peaks3['chr'] == g]
        for i in tmp3.index:
            start = tmp3.loc[i]['start']
            end = tmp3.loc[i]['end']
            mask1 = (tmp4['start'] <= end) & (tmp4['end'] >= start)
            overlap1 = tmp4[mask1]
            if len(overlap1) == 0:
                speci2 = pd.concat([speci2 , tmp3.loc[i:i]])    
                n += 1
            
    print (n)
    
    return speci2




def union_peaks(peaks1 , peaks2):
    union = peaks2.copy()
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
                pass
            else:
                union = pd.concat([union , tmp1.loc[[i]]] , axis = 0)
                
    
    union = union.drop_duplicates(subset=['chr','start','end'],keep='first')
    union = union.sort_values(by=['chr','start','end'])
    union = union.reset_index(drop=True)
    union = union[union['chr'].isin(chrom)]
    print (len(union))
    return union


def plot_venn2(n1 , n2 , n3 , title , out , lab1 , lab2):
    fig = plt.figure(figsize = (10, 10))
    venn2(subsets=(n1 , n2 , n3), set_labels=(lab1, lab2))
    for text_obj in fig.findobj(matplotlib.text.Text):
        text_obj.set_fontsize(20)
        
    plt.title(title , fontsize = 20)
    run_Plot(fig , out)
    
    
    
def plot_venn3(A_s , B_s , AB , C_s , AC , BC , ABC , title , out , lab1 , lab2 , lab3):
    fig = plt.figure(figsize = (10, 10))
    venn3(subsets=(A_s , B_s , AB , C_s , AC , BC , ABC), set_labels=(lab1, lab2 , lab3))
    for text_obj in fig.findobj(matplotlib.text.Text):
        text_obj.set_fontsize(20)
        
    plt.title(title , fontsize = 20)
    run_Plot(fig , out)
                
    
def Get_activate_SNPs(peaks , SNPs):
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

    n = 0
    act_snps = pd.DataFrame([])
    peaks.columns = ['peaks_chr' , 'peaks_start' , 'peaks_end']
    for g in chrom:
        print (g)
        tmp_snps = SNPs[SNPs['chr'] == g]
        tmp_peaks = peaks[peaks['peaks_chr'] == g]
        for i in tmp_snps.index:
            start = tmp_snps.loc[i]['start']
            end = tmp_snps.loc[i]['end']
            mask = (end >= tmp_peaks['peaks_start']) & (start <= tmp_peaks['peaks_end'])
            overlap = tmp_peaks[mask]
            if len(overlap) > 0:
                n += 1
                tmp_snps_1 = list(tmp_snps.loc[i]) + list(overlap.iloc[0])
                tmp_snps_1 = pd.DataFrame([tmp_snps_1], columns=list(tmp_snps.columns) + list(peaks.columns))
                act_snps = pd.concat([act_snps , tmp_snps_1])
    act_snps.index = range(len(act_snps))
                
    return act_snps



wt_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\peaks\\union_peaks\\HiRPC_WT_allreps_q0.05_peaks_sorted_merged.bed' , header = None)
wt_peaks.columns = ['chr' , 'start' , 'end']
ls_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\peaks\\union_peaks\\HiRPC_LS_allreps_q0.05_peaks_sorted_merged.bed' , header = None)
ls_peaks.columns = ['chr' , 'start' , 'end']
os_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\peaks\\union_peaks\\HiRPC_OS_allreps_q0.05_peaks_sorted_merged.bed' , header = None)
os_peaks.columns = ['chr' , 'start' , 'end']



#######common peaks


wt_ls_common = Common_peaks(wt_peaks , ls_peaks)
wt_os_common = Common_peaks(wt_peaks , os_peaks)
ls_os_common = Common_peaks(ls_peaks , os_peaks)



plot_venn2(len(wt_peaks) - len(wt_ls_common) , len(ls_peaks) - len(wt_ls_common) , len(wt_ls_common) , 'WT_VS_LS' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\one-dimensional\\WT_VS_LS_Hi-Coatis_peaks_overlap_Venn2.pdf' , 'WT_peaks' , 'LS_peaks')
plot_venn2(len(wt_peaks) - len(wt_os_common) , len(os_peaks) - len(wt_os_common) , len(wt_os_common) , 'WT_VS_OS' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\one-dimensional\\WT_VS_OS_Hi-Coatis_peaks_overlap_Venn2.pdf' , 'WT_peaks' , 'OS_peaks')
plot_venn2(len(ls_peaks) - len(ls_os_common) , len(os_peaks) - len(ls_os_common) , len(ls_os_common) , 'LS_VS_OS' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\one-dimensional\\LS_VS_OS_Hi-Coatis_peaks_overlap_Venn2.pdf' , 'LS_peaks' , 'OS_peaks')




common1 = wt_ls_common[[0 , 1 , 2]]
common2 = wt_os_common[[0 , 1 , 2]]
common3 = ls_os_common[[0 , 1 , 2]]

common1.columns = ['chr' , 'start' , 'end']
common2.columns = ['chr' , 'start' , 'end']
common3.columns = ['chr' , 'start' , 'end']


wt_ls_os_common = Common_peaks(common2 , ls_peaks)


wt_s = specific_peaks3(wt_peaks , ls_peaks , os_peaks)
ls_s = specific_peaks3(ls_peaks , wt_peaks , os_peaks)
os_s = specific_peaks3(os_peaks , wt_peaks , ls_peaks)






wt_ls_s = specific_peaks2(common1 , os_peaks)
wt_os_s = specific_peaks2(common2 , ls_peaks)
ls_os_s = specific_peaks2(common3 , wt_peaks)




Abc = len(wt_s)
aBc = len(ls_s)
abC = len(os_s)

ABc = len(wt_ls_s)
AbC = len(wt_os_s)
aBC = len(ls_os_s)

ABC = len(wt_ls_os_common)



out = 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\one-dimensional\\WT_VS_LS_VS_OS_Hi-Coatis_peaks_overlap_Venn3.pdf'
plot_venn3(Abc , aBc , ABc , abC , AbC , aBC , ABC , 'WT_VS_LS_VS_OS' , out , 'WT' , 'LS' , 'OS')





#####activate_SNPs

SNPs = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\CAD_related_SNPs_LD0.99_all_risk_allel_sort.narrowPeak' , header = None , usecols = (0 , 1 , 2 , 3))
SNPs.columns = ['chr' , 'start' , 'end' , 'seq']

act_snps = Get_activate_SNPs(wt_peaks , SNPs)



act_snps.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\activate_transcription_SNPs.csv' , header = True , index = None)






#####activate_regulator

nonrisk_s = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\Regulatory_Element\\HUVEC_nonrisk_silencer_q0.05_fc+-0.3_deseq2.csv' , header = 0)
nonrisk_e = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\Regulatory_Element\\HUVEC_nonrisk_enhancer_q0.05_fc+-0.3_deseq2.csv' , header = 0)


risk_s = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\Regulatory_Element\\HUVEC_risk_silencer_q0.05_fc+-0.3_deseq2.csv' , header = 0)
risk_e = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\Regulatory_Element\\HUVEC_risk_enhancer_q0.05_fc+-0.3_deseq2.csv' , header = 0)


functional_snps = list(set(list(nonrisk_s['seq']) + list(nonrisk_e['seq']) + list(risk_s['seq']) + list(risk_e['seq'])))

functional_snps_anno = []

for i in functional_snps:
    tmp1 = nonrisk_s[nonrisk_s['seq'] == i]
    tmp2 = nonrisk_e[nonrisk_e['seq'] == i]
    tmp3 = risk_s[risk_s['seq'] == i]
    tmp4 = risk_e[risk_e['seq'] == i]
    s1 = 0 ; s2 = 0 ; s3 = 0 ; s4 = 0
    if len(tmp1) > 0:
        s1 = 1
    if len(tmp2) > 0:
        s2 = 1
    if len(tmp3) > 0:
        s3 = 1
    if len(tmp4) > 0:
        s4 = 1
    
    tmp = pd.concat([tmp1 , tmp2 , tmp3 , tmp4])
    tmp_tuple = tuple(list(tmp.iloc[0]) + [s1 , s2 , s3 , s4])
    functional_snps_anno.append(tmp_tuple)
    
    
    
    
functional_snps_anno = pd.DataFrame(functional_snps_anno , columns=list(nonrisk_s.columns) + ['nonrisk_s' , 'nonrisk_e' , 'risk_s' , 'risk_e'])        






# act_nonrisk_s = Get_activate_SNPs(wt_peaks , nonrisk_s)
# act_nonrisk_e = Get_activate_SNPs(wt_peaks , nonrisk_e)



# act_risk_s = Get_activate_SNPs(wt_peaks , risk_s)
# act_risk_e = Get_activate_SNPs(wt_peaks , risk_e)



act_functional_snps_anno = Get_activate_SNPs(wt_peaks , functional_snps_anno)



speci_nonrisk_e = act_functional_snps_anno[(act_functional_snps_anno['nonrisk_e'] == 1) & (act_functional_snps_anno['risk_e'] == 0)]
speci_nonrisk_s = act_functional_snps_anno[(act_functional_snps_anno['nonrisk_s'] == 1) & (act_functional_snps_anno['risk_s'] == 0)]


speci_risk_e = act_functional_snps_anno[(act_functional_snps_anno['nonrisk_e'] == 0) & (act_functional_snps_anno['risk_e'] == 1)]
speci_risk_s = act_functional_snps_anno[(act_functional_snps_anno['nonrisk_s'] == 0) & (act_functional_snps_anno['risk_s'] == 1)]


speci_nonrisk_e.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\activate_regulatory_elements\\HUVEC_activate_nonrisk_specific_enhancer.csv' , header = True , index = None)
speci_nonrisk_s.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\activate_regulatory_elements\\HUVEC_activate_nonrisk_specific_silencer.csv' , header = True , index = None)
speci_risk_e.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\activate_regulatory_elements\\HUVEC_activate_risk_specific_enhancer.csv' , header = True , index = None)
speci_risk_s.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\activate_regulatory_elements\\HUVEC_activate_risk_specific_silencer.csv' , header = True , index = None)








##############Hi-Coatis_VS_ATAC
wt_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_WT_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
wt_atac_peak.columns = ['chr' , 'start' , 'end']
ls_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_LS_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
ls_atac_peak.columns = ['chr' , 'start' , 'end']
os_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_OS_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
os_atac_peak.columns = ['chr' , 'start' , 'end']



wt_h_a_common = Common_peaks(wt_peaks , wt_atac_peak)
ls_h_a_common = Common_peaks(ls_peaks , ls_atac_peak)
os_h_a_common = Common_peaks(os_peaks , os_atac_peak)



plot_venn2(len(wt_peaks) - len(wt_h_a_common) , len(wt_atac_peak) - len(wt_h_a_common) , len(wt_h_a_common) , 'WT_HiCoatis_VS_ATAC' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig2\\WT_HiCoatis_VS_ATAC_peaks_overlap_Venn2.pdf' , 'WT_HiCoatis' , 'WT_ATAC')
plot_venn2(len(ls_peaks) - len(ls_h_a_common) , len(ls_atac_peak) - len(ls_h_a_common) , len(ls_h_a_common) , 'LS_HiCoatis_VS_ATAC' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig2\\LS_HiCoatis_VS_ATAC_peaks_overlap_Venn2.pdf' , 'LS_HiCoatis' , 'LS_ATAC')
plot_venn2(len(os_peaks) - len(os_h_a_common) , len(os_atac_peak) - len(os_h_a_common) , len(os_h_a_common) , 'OS_HiCoatis_VS_ATAC' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig2\\OS_HiCoatis_VS_ATAC_peaks_overlap_Venn2.pdf' , 'OS_HiCoatis' , 'OS_ATAC')





##############specific Hi-Coatis_VS_ATAC
speci_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_specific_peaks_edgR_q0.05_fc1.bed' , header = None)
speci_ls.columns = ['chr' , 'start' , 'end']
speci_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_OS_specific_peaks_edgR_q0.05_fc1.bed' , header = None)
speci_os.columns = ['chr' , 'start' , 'end']
common = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_OS_common_peaks_edgR_q0.05_fc1.bed' , header = None)
common.columns = ['chr' , 'start' , 'end']


wt_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_WT_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
wt_atac_peak.columns = ['chr' , 'start' , 'end']
ls_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_LS_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
ls_atac_peak.columns = ['chr' , 'start' , 'end']
os_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_OS_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
os_atac_peak.columns = ['chr' , 'start' , 'end']
ls_os_union_atac_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_LS_OS_union_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
ls_os_union_atac_peaks.columns = ['chr' , 'start' , 'end']



speci_ls_common = Common_peaks(speci_ls , ls_atac_peak)
speci_os_common = Common_peaks(speci_os , os_atac_peak)
common_common = Common_peaks(common , ls_os_union_atac_peaks)



plot_venn2(len(speci_ls) - len(speci_ls_common) , len(ls_atac_peak) - len(speci_ls_common) , len(speci_ls_common) , 'LS_speci_HiCoatis_VS_ATAC' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig3\\LS_speci_HiCoatis_VS_ATAC_peaks_overlap_Venn2.pdf' , 'LS_speci_HiCoatis' , 'LS_ATAC')
plot_venn2(len(speci_os) - len(speci_os_common) , len(os_atac_peak) - len(speci_os_common) , len(speci_os_common) , 'OS_speci_HiCoatis_VS_ATAC' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig3\\OS_speci_HiCoatis_VS_ATAC_peaks_overlap_Venn2.pdf' , 'OS_speci_HiCoatis' , 'OS_ATAC')





speci_ls_common[[0 , 1 , 2]].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\specific_peaks_Anno\\HUVEC_LS_specific_peaks_edgR_q0.05_fc1_overlaped_ATAC.bed' , header = None , index = None , sep = '\t')
speci_os_common[[0 , 1 , 2]].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\specific_peaks_Anno\\HUVEC_OS_specific_peaks_edgR_q0.05_fc1_overlaped_ATAC.bed' , header = None , index = None , sep = '\t')







