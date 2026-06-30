# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 17:07:11 2025

@author: lenovo
"""


import numpy as np
import pandas as pd




chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

def pairs_to_ipet(file , out , headline):
    data = pd.read_csv(file , header = None , usecols=(1 , 2 , 3, 4 , 5 , 6), skiprows = headline , sep = '\t')
    
    data.columns = ['chr1' , 'pos1' , 'chr2' , 'pos2' , 'strand1' , 'strand2']
    ##len(data) = 374649585
    
    
    data_new = data[['chr1' , 'pos1' , 'strand1' , 'chr2' , 'pos2' , 'strand2']]
    
    
    
    data_new = data_new[(data_new['pos2'] - data_new['pos1']) >= 8000]
    ##len(data_new) = 73731410
    data_new = data_new[data_new['chr1'] == data_new['chr2']]
    ##len(data_new) = 50574413
    data_new = data_new[data_new['chr1'].isin(chrom)]
    data_new = data_new[data_new['chr2'].isin(chrom)]
    ##len(data_new) = 50328225
    
    data_new.to_csv(out , header=None , index = None , sep = '\t')
    
    
    
    



    
    
#############HUVEC_control_LS_OS
 
pairs_to_ipet('/scratch/2025-12-01/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/Combined_HUVEC_control_merged4.hg38.nodups.pairs_+8kb' , \
              '/scratch/2025-12-01/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/loops/HUVEC_WT/ChIA-pipline_8K+/Combined_HUVEC_control_merged4.hg38.nodups.pairs_+8kb.ipet' , \
              943) 
    
    
pairs_to_ipet('/scratch/2025-12-01/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/Combined_HUVEC_ls_merged4.hg38.nodups.pairs_+8kb' , \
              '/scratch/2025-12-01/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/loops/HUVEC_LS/ChIA-pipline_8K+/Combined_HUVEC_ls_merged4.hg38.nodups.pairs_+8kb.ipet' , \
              943) 
    
    
    
pairs_to_ipet('/scratch/2025-12-01/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/Combined_HUVEC_os_merged8.hg38.nodups.pairs_+8kb' , \
              '/scratch/2025-12-01/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/loops/HUVEC_OS/ChIA-pipline_8K+/Combined_HUVEC_os_merged8.hg38.nodups.pairs_+8kb.ipet' , \
              971) 
    
    

    
##############WZ


    
pairs_to_ipet('/scratch/2025-12-08/bio-shenw/WZ/Hi_RPC/workspace_CRISPR/results/pairs_library/KHM2_Rep1.hg38.nodups.pairs_+8kb' , \
              '/scratch/2025-12-08/bio-shenw/WZ/Hi_RPC/workspace_CRISPR/results/loops/KHM2/ChIA-pipline_8K+/KHM2_Rep1.hg38.nodups.pairs_+8kb.ipet' , \
              921) 
    
    
    
    
    
    
    
#############HUVEC_OS_allreps
 
pairs_to_ipet('/scratch/2026-06-15/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC_allReps/Combined_huvec_os_merged11.hg38.nodups.pairs_+8kb' , \
              '/scratch/2026-06-15/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC_allReps/Combined_huvec_os_merged11.hg38.nodups.pairs_+8kb.ipet' , \
              1003) 
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

