# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 21:24:20 2026

@author: lenovo
"""

from __future__ import division
from scipy import sparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import sys
import os
from skimage.transform import resize



# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

# Tad_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
#                     'formats':['U8' , np.int , np.int]})

def normalize_TAD_interaction(Lib , Tad , length = 40):
    
    Tad = Tad[(Tad['end'] - Tad['start']) >= R * 10]
    IS = np.zeros(length + 20)
    resize_matrix = np.zeros((length , length))
    n = 0
    for i in Tad.index:
        chro = Tad.loc[i]['chr']
        start = Tad.loc[i]['start'] // R 
        end = Tad.loc[i]['end'] // R
        d = end - start 
        # if end - start < 20:
        #     continue
        startHiC = start - d 
        endHiC = end + d 
        if startHiC < 0:
            continue
        matrix = Lib[chro][startHiC:endHiC , startHiC:endHiC]
        resize_m = resize(matrix , (length + 20 , length + 20), order = 3)
        resize_matrix += resize_m[10:50,10:50] 
        n += 1
        for bound in range(length + 20):
            Is = acquireSingleIns(resize_m , bound  , 10 , 'divide')
            IS[bound] += Is
            
    print (n)
    resize_matrix = resize_matrix / n
    IS = IS / n
    return resize_matrix, IS[10:50]
    
    
    
def acquireSingleIns(matrix_data_chr,bound,left_right,category):
    ins=0
    start_site=0;end_site=matrix_data_chr.shape[0]
    if ((bound-left_right<start_site)|(end_site<bound+left_right)):        
        return ins    
    aa=matrix_data_chr[bound-left_right:bound-0,bound+0:bound+left_right]
    b1=[[matrix_data_chr[i,j] for i in range(bound-left_right,bound-0) if j>i] 
            for j in range(bound-left_right,bound-0)]
    b2=[[matrix_data_chr[i,j] for i in range(bound+0,bound+left_right) if j>i] 
            for j in range(bound+0,bound+left_right)]
    
    aa_zero=sum([sum(np.array(item)==0) for item in aa])
    b1_zero=sum([sum(np.array(item)==0) for item in b1])
    b2_zero=sum([sum(np.array(item)==0) for item in b2])
    if aa_zero+b1_zero+b2_zero>=left_right:
        return ins    
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
    
def Heatmap_Plot(matrix ,  c):
    matrix = matrix / matrix.sum()
    vmax = np.percentile(matrix , 85)
    vmin = np.percentile(matrix , 3)
    
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                    extent = (0, len(matrix), 0, len(matrix)), vmax = 0.0009411, vmin = 0.0000728 , origin = 'lower')
    cxlim = ax1.get_xlim()
    cylim = ax1.get_ylim()
    ## Ticks and Labels
    ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
    labels = ['-1/2TAD' , 'Boundary' , '' , 'Boundary' , '+1/2TAD']
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)
    ax1.set_yticks(ticks)
    ax1.set_yticklabels(labels, rotation = 'horizontal')
    ax1.set_xlabel(c , fontsize = 30 )
    
    ax1.set_xlim(cxlim)
    ax1.set_ylim(cylim)                    
    ## Colorbar
    ax2 = fig.add_axes([Left + 0.6 , HB - 0.06 , 0.1 , 0.035])
    fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [0.0000728 ,0.0009411])
    return fig

def Insulation_Plot(data):
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    ax1.plot(np.arange(len(data[0])) , data[0] , label = 'ST')
    ax1.plot(np.arange(len(data[1])) , data[1] , label = 'LSS')
    ax1.plot(np.arange(len(data[2])) , data[2] , label = 'OSS')

    
    ## Ticks and Labels
#    ticks = list(np.linspace(0 , len(matrix) , 7).astype(float))
#    labels = ['-TAD' , '-1/2TAD' , 'Boundary' , '' , 'Boundary' , '+1/2TAD' , '+TAD']
#    ax1.set_xticks(ticks)
#    ax1.set_xticklabels(labels)
#    ax1.set_yticks(ticks)
#    ax1.set_yticklabels(labels, rotation = 'horizontal')
    ax1.set_ylabel('Insulation Score' , fontsize = 20 )
    ax1.legend()
    

    return fig

    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
        
    
def Matrix_nan_to_0(data):
    Matrix = {}
    for g in data:
        matrix = data[g]
        matrix[np.isnan(matrix)] = 0
        Matrix[g] = matrix
    return Matrix
    

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

R = 25000
size = (11, 12)
Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7



##Load HiC_data
wt = np.load('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\npz\\HUVEC_WT_matrix_diagonal_0_dict_25K.npz')
lss = np.load('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\npz\\HUVEC_LSS_matrix_diagonal_0_dict_25K.npz')
oss = np.load('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\npz\\HUVEC_OSS_matrix_diagonal_0_dict_25K.npz')


    
    


##Load Tad_data
# CCS_Tad = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\CCS_Domain_bottom_40K_respective_stable_400K.txt' , dtype = Tad_type)
# CCS_Tad = CCS_Tad[CCS_Tad['chr'] != 'X']
# Tad = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\union_domain_NTs_fESC_merged.txt' , dtype = Tad_type)
# Tad = Tad[Tad['chr'] != 'X']

Tad = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\union_TADs_25K_domains.bed' , header = None)
Tad.columns = ['chr' , 'start' , 'end']




matrix_wt , IS_wt = normalize_TAD_interaction(wt , Tad , 40)
matrix_lss , IS_lss = normalize_TAD_interaction(lss , Tad , 40)
matrix_oss , IS_oss = normalize_TAD_interaction(oss , Tad , 40)



fig_wt = Heatmap_Plot(matrix_wt , 'ST')
fig_lss = Heatmap_Plot(matrix_lss , 'LSS')
fig_oss = Heatmap_Plot(matrix_oss , 'OSS')




run_Plot(fig_wt , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\Tad_normalization_interaction_plot\\WT_average_Tad_interaction_25K.pdf')
run_Plot(fig_lss , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\Tad_normalization_interaction_plot\\LSS_average_Tad_interaction_25K.pdf')
run_Plot(fig_oss , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\Tad_normalization_interaction_plot\\OSS_average_Tad_interaction_25K.pdf')


fig_IS = Insulation_Plot([IS_wt , IS_lss , IS_oss])
run_Plot(fig_IS , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\Tad_normalization_interaction_plot\\Tad_Insulation_score_25K.pdf')




                    