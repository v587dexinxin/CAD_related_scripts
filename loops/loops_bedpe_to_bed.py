# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 16:14:27 2026

@author: lenovo
"""

import matplotlib.pyplot as plt
import pandas as pd




loops_wt = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bedpe' , header = None , usecols = (0 , 1 , 2 , 3 , 4 , 5))
loops_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None , usecols = (0 , 1 , 2 , 3 , 4 , 5))
loops_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None , usecols = (0 , 1 , 2 , 3 , 4 , 5))

loops_wt.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2']
loops_ls.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2']
loops_os.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2']





def Get_center_loops(loops):
    loops['START'] = (loops['end1'] + loops['start1']) // 2
    loops['END'] = (loops['end2'] + loops['start2']) // 2
    return loops
    



loops_wt = Get_center_loops(loops_wt)
loops_ls = Get_center_loops(loops_ls)
loops_os = Get_center_loops(loops_os)



loops_wt[['chr1' , 'START' , 'END']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bed' , header = None , index = None , sep = '\t')
loops_ls[['chr1' , 'START' , 'END']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bed' , header = None , index = None , sep = '\t')
loops_os[['chr1' , 'START' , 'END']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bed' , header = None , index = None , sep = '\t')






















