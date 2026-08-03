# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 15:32:49 2026

@author: lenovo
"""

import pandas as pd
import numpy as np



#########loops#############
loops_wt = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops_wt.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']
loops_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops_ls.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']
loops_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops_os.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']




#####selected_loops_s42
s42_loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\OS_s42_chr11_118910720_118913366_peaks_+-1kb_seq2178_rs11606719_118911765_HUVEC_nonrisk_enhancer_speci.bedpe' , header = None)
s42_loops.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']
RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)

RNA['Start'] = RNA['Start'] - 2000
RNA['End'] = RNA['End'] + 2000



rna_new = pd.DataFrame([])

n = 1
for i in s42_loops.index:
    print ('loop' + str(n))
    g = s42_loops.iloc[i]['chr1']
    s1 = s42_loops.iloc[i]['s1']
    e1 = s42_loops.iloc[i]['e1']
    s2 = s42_loops.iloc[i]['s2']
    e2 = s42_loops.iloc[i]['e2']
    tmp_rna = RNA[RNA['Chr'] == g]
    mask1 = (tmp_rna['Start'] <= e1) & (tmp_rna['End'] >= s1)
    mask2 = (tmp_rna['Start'] <= e2) & (tmp_rna['End'] >= s2)
    overlap1 = tmp_rna[mask1]
    overlap2 = tmp_rna[mask2]
    if len(overlap1) > 0:
        rna_new = pd.concat([rna_new , overlap1])
    else:
        pass
    if len(overlap2) > 0:
        rna_new = pd.concat([rna_new , overlap2])
    else:
        pass
    n += 1
    
    
rna_new = rna_new.drop_duplicates(rna_new)    




#####selected_loops_seq9454
s42_loops = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_loops\\WT_seq9454.bedpe' , header = None)
s42_loops.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']
RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)

RNA['Start'] = RNA['Start'] - 2000
RNA['End'] = RNA['End'] + 2000



rna_new = pd.DataFrame([])

n = 1
for i in s42_loops.index:
    print ('loop' + str(n))
    g = s42_loops.iloc[i]['chr1']
    s1 = s42_loops.iloc[i]['s1']
    e1 = s42_loops.iloc[i]['e1']
    s2 = s42_loops.iloc[i]['s2']
    e2 = s42_loops.iloc[i]['e2']
    tmp_rna = RNA[RNA['Chr'] == g]
    mask1 = (tmp_rna['Start'] <= e1) & (tmp_rna['End'] >= s1)
    mask2 = (tmp_rna['Start'] <= e2) & (tmp_rna['End'] >= s2)
    overlap1 = tmp_rna[mask1]
    overlap2 = tmp_rna[mask2]
    if len(overlap1) > 0:
        rna_new = pd.concat([rna_new , overlap1])
    else:
        pass
    if len(overlap2) > 0:
        rna_new = pd.concat([rna_new , overlap2])
    else:
        pass
    n += 1
    
    
rna_new = rna_new.drop_duplicates(rna_new)    












