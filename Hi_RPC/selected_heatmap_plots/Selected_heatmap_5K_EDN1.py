# -*- coding: utf-8 -*-
"""
Created on Thu May 21 17:37:35 2026

@author: lenovo
"""

from __future__ import division
import numpy as np
import cooler
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF' ,'#CD0000'])
my_cmap.set_bad('#D3D3D3')




def nan_to_zero(matrix):
    '''
    '''
    nanmask = np.isnan(matrix)
    matrix[nanmask] = 0
    return matrix

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'M'])
    else:
        return ''.join([str(i_part), 'M'])
    

def Get_matrix(file_name , g , res):
    '''
    '''
    c = cooler.Cooler(file_name + "::/resolutions/" + res)
    c_matrix = c.matrix(balance=False).fetch(g)
    matrix = nan_to_zero(c_matrix)
    
    return matrix
    
    

####matrix#######


f_wt = "/scratch/2026-06-29/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_WT_merged.hg38.nodups.mapq_30.1000.mcool"
f_ls = "/scratch/2026-06-29/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_LS_merged.hg38.nodups.mapq_30.1000.mcool"
f_os = "/scratch/2026-06-29/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_OS_merged.hg38.nodups.mapq_30.1000.mcool"

g = 'chr6'
R = 5000
gene_name = 'EDN1'

interval = [(g , 11300000 , 12800000)]


startHiC = interval[0][1] // R 
endHiC = interval[0][2] // R 



wt_matrix = Get_matrix(f_wt, g, str(R))[startHiC:endHiC , startHiC:endHiC]
ls_matrix = Get_matrix(f_ls, g, str(R))[startHiC:endHiC , startHiC:endHiC]
os_matrix = Get_matrix(f_os, g, str(R))[startHiC:endHiC , startHiC:endHiC]



data = {'WT' : wt_matrix , 'LSS' : ls_matrix , 'OSS' : os_matrix}

##

cells = ['WT' , 'LSS' , 'OSS']

for c in cells:
    print (c)
    
    matrix_0 = data[c]
    
    np.fill_diagonal(matrix_0, 0)
    
    
    #####Out_Files#####
    OutFolder = '/scratch/2026-06-29/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/selected_heatmap_plots_' + gene_name
    OutFil = 'HUVEC_' + c + '_heatmap_' + g + '_' + str(R // 1000) + 'K_' + gene_name + '.pdf'
    pp = PdfPages(os.path.join(OutFolder , OutFil))
    
    

    
    
    ####Plot####
    size = (12, 12)   
    Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
    

    
    matrix = matrix_0
    
    ticks = list(np.linspace(0 , matrix.shape[1] , 2).astype(float))
    pos = [((startHiC + t) * R) for t in ticks]
    labels = [properU(p) for p in pos]
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero, 55)
    if c == 'LSS':
        vmax = np.percentile(nonzero, 50)
    elif c == 'OSS':
        vmax = np.percentile(nonzero, 10)
    else:
        pass

    
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left  , HB , width , HH])
    sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none', vmax = vmax,extent = (0, matrix.shape[1], matrix.shape[0] , 0) , origin = 'upper')
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(labels , fontsize=15)
    ax.set_yticklabels(labels , fontsize=15)
    ax.set_xlabel(c + '_' + g , fontsize=20)
    ax.yaxis.set_label_position("right")
    # ax.set_xlim(0, matrix.shape[1])
    # ax.set_ylim(0, matrix.shape[0])
    ax.set_title('HUVEC_' + c + '_heatmap_' + g + '_' + str(R // 1000) + 'K' , fontsize=25)
    
    # ##Selected_interval
    
    
    
    ## Colorbar
    ax = fig.add_axes([Left + 0.5 , HB - 0.1 , 0.1 , 0.035])
    cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
    cbar.set_ticks([0 , vmax])
    
    
    
    
    pp.savefig(fig)
    pp.close()