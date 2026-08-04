# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:41:01 2020

@author: xxli
"""


from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import sys
import os

import matplotlib.pyplot as plt



 
def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in HiC_Lib:
        tmp = HiC_Lib[g]
        tmp[np.isnan(tmp)] = 0
        Lib_new[g] = tmp
    return Lib_new
                               

def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = True, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , length = 5 ,labelsize = 30  )

def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
def getmatrix(inter,l_bin,r_bin):
    inter_matrix = np.zeros((r_bin - l_bin, r_bin - l_bin),dtype = float )
    mask = (inter['bin1'] >= l_bin) & (inter['bin1'] < r_bin) & \
           (inter['bin2'] >= l_bin) & (inter['bin2'] < r_bin)
    inter_extract = inter[mask]
    for i in inter_extract:
        if i['bin1'] != i['bin2']:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
            inter_matrix[i['bin2'] - l_bin][i['bin1'] - l_bin] += i['IF']
        else:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
    return inter_matrix

def OE_matrix(matrix):
    matrix_new = np.zeros((len(matrix) , len(matrix)))
    exp = []
    for i in range(len(matrix)):
        tmp = []
        for j in range(len(matrix) - i):
            tmp.append((matrix[j][j + i]))
        m = sum(tmp)/len(tmp)
        if m == 0:
            exp.append(1)
        else:
            exp.append(m)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            dis = abs(i-j)
            if matrix[i][j] == 0:
                matrix_new[i][j] = 0
            else:
                matrix_new[i][j] = np.log2(matrix[i][j] / exp[dis])
    return matrix_new
        

def acquireSingleIns(matrix_data_chr,bound,left_right,category):
    ins=0
    start_site=0;end_site=matrix_data_chr.shape[0]
    if ((bound-left_right<start_site)|(end_site<bound+left_right)):        
        return ins    

    aa=matrix_data_chr[bound-left_right+1:bound+1,bound+1:bound+left_right+1]
    b1=[[matrix_data_chr[i,j] for i in range(bound-left_right,bound) if j>i] 
            for j in range(bound-left_right,bound)]
    b2=[[matrix_data_chr[i,j] for i in range(bound+1,bound+left_right+1) if j>i] 
            for j in range(bound+1,bound+left_right+1)]
    
    aa_zero=sum([sum(np.array(item)==0) for item in aa])
    b1_zero=sum([sum(np.array(item)==0) for item in b1])
    b2_zero=sum([sum(np.array(item)==0) for item in b2])
    # if aa_zero+b1_zero+b2_zero>=left_right:
    #     return ins    
    aa_sum=sum([sum(item) for item in aa])
    b1_sum=sum([sum(item) for item in b1])
    b2_sum=sum([sum(item) for item in b2])
    if aa_sum>0:
        if(category=='divide'):
            ins=np.log2((aa_sum+b1_sum+b2_sum)/float(aa_sum))
        elif(category=='average'):
            ins=aa_sum/float(left_right)/left_right
        else:
            print('the calc type went wrong')
    return ins

    
                
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['U8' , np.int , np.int]})

boundary_type = np.dtype({'names':['chr' , 'pos'],
                     'formats':['U8' , np.int]})


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
cells = ['ST' , 'LSS' , 'OSS']



#HiC Data Process
WT_Lib = np.load('/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/npz/HUVEC_WT_matrix_diagonal_0_dict_25K.npz')
WT_Lib = Get_nan_zero_Matrix(WT_Lib)
LSS_Lib = np.load('/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/npz/HUVEC_LSS_matrix_diagonal_0_dict_25K.npz')
LSS_Lib = Get_nan_zero_Matrix(LSS_Lib)
OSS_Lib = np.load('/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/npz/HUVEC_OSS_matrix_diagonal_0_dict_25K.npz')
OSS_Lib = Get_nan_zero_Matrix(OSS_Lib)



HiC_Data = {'ST' : WT_Lib , 'LSS' : LSS_Lib , 'OSS' : OSS_Lib}

 

#-------------------------------------insulation_score-----------------------------------------           
Tad = np.loadtxt('/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/h5/TADs/res_10K/union_boundary_10K.bed' , dtype = boundary_type)
#Tad = Tad[Tad['chr'] != 'X']


for c in cells:
    out = open('/scratch/2026-07-06/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/h5/TADs/res_10K/Insulation_score/' + c + '_Insulation_score_25K_15.txt' , 'w')
    out.writelines('\t'.join(['chr' , 'Boundary' , 'Insulation_Score']) + '\n')
    HiC_Lib = HiC_Data[c]
    Tad_Lib = Tad
    for i in Tad_Lib:
        chro = i['chr']
        boundary = i['pos']
        startHiC = (boundary - 500000) // 25000
        endHiC = (boundary + 500000) // 25000
        HiC_lib = HiC_Lib[chro]
        if endHiC > len(HiC_lib):
            continue
        matrix = HiC_lib[startHiC:endHiC , startHiC:endHiC]
        insulation = acquireSingleIns(matrix , 20 , 15 , 'divide')
        out.writelines('\t'.join([chro , str(boundary) , str(insulation)]) + '\n')
        
    out.close()
        


    