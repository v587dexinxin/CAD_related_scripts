# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 14:30:33 2025

@author: lenovo
"""


from __future__ import division
import numpy as np 
import pandas as pd



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
            if overlap.size != 0:
                n += 1
                common.append((g , start , end))
    print (n)
    common = pd.DataFrame(common)
    common.columns = ['chr' , 'start' , 'end']
    
    return common










##################HUVEC_WT

loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Combined_HUVEC_control_merged4.hg38.nodups.pairs_+8kb.cluster.FDRfiltered.txt' , usecols = (0 , 1 , 2 , 3 , 4 , 5 , 6 , 11) , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'pvalue']


peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_WT_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks.columns = ['chr' , 'start' , 'end']


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


common = pd.DataFrame([])

for g in chrom:
    print(g)
    tmp_peaks = peaks[peaks['chr'] == g]
    tmp_loops = loops[loops['chr1'] == g]
    for i in tmp_loops.index:
        start1 = tmp_loops.loc[i]['start1']
        end1 = tmp_loops.loc[i]['end1']
        start2 = tmp_loops.loc[i]['start2']
        end2 = tmp_loops.loc[i]['end2']
        mask1 = (tmp_peaks['start'] <= end1) & (tmp_peaks['end'] >= start1)
        mask2 = (tmp_peaks['start'] <= end2) & (tmp_peaks['end'] >= start2)
        overlap1 = tmp_peaks[mask1]
        overlap2 = tmp_peaks[mask2]
        if abs(start2 - start1) >= 0:
            if (overlap1.size != 0) or (overlap2.size != 0):
                common = pd.concat([common , tmp_loops.loc[i:i]] , axis = 0)
        

common = common.drop_duplicates()
    
common.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bedpe' , header = None , index = None , sep = '\t')





##################HUVEC_LS

loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Combined_HUVEC_ls_merged4.hg38.nodups.pairs_+8kb.cluster.FDRfiltered.txt' , usecols = (0 , 1 , 2 , 3 , 4 , 5 , 6 , 11) , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'pvalue']


peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_LS_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks.columns = ['chr' , 'start' , 'end']


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


common = pd.DataFrame([])

for g in chrom:
    print(g)
    tmp_peaks = peaks[peaks['chr'] == g]
    tmp_loops = loops[loops['chr1'] == g]
    for i in tmp_loops.index:
        start1 = tmp_loops.loc[i]['start1']
        end1 = tmp_loops.loc[i]['end1']
        start2 = tmp_loops.loc[i]['start2']
        end2 = tmp_loops.loc[i]['end2']
        mask1 = (tmp_peaks['start'] <= end1) & (tmp_peaks['end'] >= start1)
        mask2 = (tmp_peaks['start'] <= end2) & (tmp_peaks['end'] >= start2)
        overlap1 = tmp_peaks[mask1]
        overlap2 = tmp_peaks[mask2]
        if abs(start2 - start1) >= 0:
            if (overlap1.size != 0) or (overlap2.size != 0):
                common = pd.concat([common , tmp_loops.loc[i:i]] , axis = 0)
        

common = common.drop_duplicates()
    
common.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None , index = None , sep = '\t')




##################HUVEC_OS

loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Combined_HUVEC_os_merged8.hg38.nodups.pairs_+8kb.cluster.FDRfiltered.txt' , usecols = (0 , 1 , 2 , 3 , 4 , 5 , 6 , 11) , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'pvalue']


peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_OS_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks.columns = ['chr' , 'start' , 'end']


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


common = pd.DataFrame([])

for g in chrom:
    print(g)
    tmp_peaks = peaks[peaks['chr'] == g]
    tmp_loops = loops[loops['chr1'] == g]
    for i in tmp_loops.index:
        start1 = tmp_loops.loc[i]['start1']
        end1 = tmp_loops.loc[i]['end1']
        start2 = tmp_loops.loc[i]['start2']
        end2 = tmp_loops.loc[i]['end2']
        mask1 = (tmp_peaks['start'] <= end1) & (tmp_peaks['end'] >= start1)
        mask2 = (tmp_peaks['start'] <= end2) & (tmp_peaks['end'] >= start2)
        overlap1 = tmp_peaks[mask1]
        overlap2 = tmp_peaks[mask2]
        if abs(start2 - start1) >= 0:
            if (overlap1.size != 0) or (overlap2.size != 0):
                common = pd.concat([common , tmp_loops.loc[i:i]] , axis = 0)
        

common = common.drop_duplicates()
    
common.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops_1.bedpe' , header = None , index = None , sep = '\t')








##################HUVEC_OS_chr4

loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Combined_HUVEC_os_merged8.hg38.nodups.pairs_+8kb.cluster.FDRfiltered.txt' , usecols = (0 , 1 , 2 , 3 , 4 , 5 , 6 , 11) , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'pvalue']


peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_OS_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks.columns = ['chr' , 'start' , 'end']


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


common = pd.DataFrame([])

for g in ['chr4']:
    print(g)
    tmp_peaks = peaks[peaks['chr'] == g]
    tmp_loops = loops[loops['chr1'] == g]
    for i in tmp_loops.index:
        start1 = tmp_loops.loc[i]['start1']
        end1 = tmp_loops.loc[i]['end1']
        start2 = tmp_loops.loc[i]['start2']
        end2 = tmp_loops.loc[i]['end2']
        mask1 = (tmp_peaks['start'] <= end1) & (tmp_peaks['end'] >= start1)
        mask2 = (tmp_peaks['start'] <= end2) & (tmp_peaks['end'] >= start2)
        overlap1 = tmp_peaks[mask1]
        overlap2 = tmp_peaks[mask2]
        if abs(start2 - start1) >= 0:
            if (overlap1.size != 0) or (overlap2.size != 0):
                common = pd.concat([common , tmp_loops.loc[i:i]] , axis = 0)
        

common = common.drop_duplicates()
    
common.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops_chr4.bedpe' , header = None , index = None , sep = '\t')







