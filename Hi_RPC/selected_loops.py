# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 20:18:56 2025

@author: lenovo
"""

import pandas as pd
import numpy as np
from itertools import islice

loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' ,'IF' , 'pvalue']

tmp = loops[(loops['chr1'] == 'chr9') & (loops['chr2'] == 'chr9')]

'''
DIRC3: chr2: 217282019-217790443
'''

start = 107000000
end = 109000000


dirc3 = tmp[((tmp['s1'] >= start) & (tmp['e1'] <= end)) & ((tmp['s2'] >= start) & (tmp['e2'] <= end))]

dirc3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Selected_Klf4_HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None , index = None , sep = '\t')





loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Selected_Klf4_HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
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
            if (overlap1.size != 0) and (overlap2.size != 0):
                common = pd.concat([common , tmp_loops.loc[i:i]] , axis = 0)
        

common = common.drop_duplicates()

common.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\Selected_Klf4_HUVEC_LS_HiCoatis_twoanchor_binding_loops.bedpe' , header = None , index = None , sep = '\t')

