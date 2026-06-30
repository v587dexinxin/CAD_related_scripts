# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 22:28:45 2024

@author: lenovo
"""


from __future__ import division
import numpy as np
import pandas as pd
#from tadlib.calfea.analyze import getmatrix
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys



def Box_plot_3cellline(data , vmin , vmax , title):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'green','linewidth':2},
            medianprops={'color':'green','linewidth':2},
            capprops={'color':'green','linewidth':2},
            whiskerprops={'color':'green','linewidth':2})
    
    for i, value in enumerate(data[0]):
        ax.plot(1, value, 'ro' , c = 'darkred')  # 在第一个箱线图上标注数据点
    
    for i, value in enumerate(data[1]):
        ax.plot(2, value, 'ro' , c = 'dodgerblue')  # 在第一个箱线图上标注数据点
        
    for i, value in enumerate(data[2]):
        ax.plot(3, value, 'ro' , c = 'green')  # 在第一个箱线图上标注数据点
    
    
    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)
    
    
    ax.set_xticks([1 , 2 , 3 ])
    ax.set_xticklabels(['HUVEC_ST' , 'HUVEC_LS' , 'HUVEC_OS'] , fontsize = 20)
    ax.set_ylabel('Gene Expression (log2(FPKM + 1)' , fontsize = 20)
    ax.set_xlabel(title , fontsize = 20)
    ax.set_xlim((0.5 , 3.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()





# data = pd.read_excel('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA_seq_Luhaocheng\\副本genecode_v23_STvsLS_ge4with1fpm_limma.xlsx' , header = 0) 



union = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0 , index_col=0)
union['ST_FPKM'] = (union['ST_R1_FPKM'] + union['ST_R2_FPKM']) / 2
union['LS_FPKM'] = (union['LS_R1_FPKM'] + union['LS_R2_FPKM']) / 2
union['OS_FPKM'] = (union['OS_R1_FPKM'] + union['OS_R2_FPKM']) / 2
data = union.copy()





TRIB1 = data.loc['TRIB1']

x = list(TRIB1[['ST_R1_FPKM', 'ST_R2_FPKM']])

y = list(TRIB1[['LS_R1_FPKM', 'LS_R2_FPKM']])

z = list(TRIB1[['OS_R1_FPKM', 'OS_R2_FPKM']])

fig = Box_plot_3cellline([x , y , z] , 0 , 6 , 'TRIB1 gene expression')


run_Plot(fig , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\HUVEC_control_LS_TRIB1_geneexpression_boxplot_new.pdf')