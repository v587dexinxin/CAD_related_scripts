# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 21:29:23 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import cooler






def nan_to_zero(matrix):
    '''
    '''
    nanmask = np.isnan(matrix)
    matrix[nanmask] = 0
    return matrix



def Get_matrix(file_name , res):
    '''
    '''
    k = 2
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    matrix_new = {}
    for g in chrom:
        c = cooler.Cooler(file_name + "::/resolutions/" + res)
        c_matrix = c.matrix(balance=True).fetch(g)
        matrix = nan_to_zero(c_matrix)
        
        ##########对角线为0#########
        for i in range(-k, k + 1):
            rows, cols = np.diag_indices_from(matrix)
        
            if i > 0:
                matrix[rows[:-i], cols[i:]] = 0
            elif i < 0:
                matrix[rows[-i:], cols[:i]] = 0
            else:
                matrix[rows, cols] = 0
                
        matrix_new[g] = matrix
    return matrix_new




f_wt = "/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_WT_merged.hg38.nodups.mapq_30.1000.mcool"
f_ls = "/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_LS_merged.hg38.nodups.mapq_30.1000.mcool"
f_os = "/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_OS_merged.hg38.nodups.mapq_30.1000.mcool"


R = 25000
    
    
wt_matrix = Get_matrix(f_wt, str(R))
ls_matrix = Get_matrix(f_ls, str(R))
os_matrix = Get_matrix(f_os, str(R))


data = {'WT' : wt_matrix , 'LSS' : ls_matrix , 'OSS' : os_matrix}





for c in ['WT' , 'LSS' , 'OSS']:
    np.savez('HUVEC_' + c + '_matrix_diagonal_0_dict_' + str(R // 1000) + 'K.npz' , **data[c])

















































