# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:11:03 2024

@author: lenovo
"""


from __future__ import division
import numpy as np 
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
import scipy
from scipy.stats import ranksums
from matplotlib_venn import venn2, venn2_circles
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

def Box_plot_4cellline(data , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=True, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=True, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[2] , positions=[4] , showfliers=True, widths = 0.7 ,
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[3] , positions=[5] , showfliers=True, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[2] , data[3])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[2] , data[3])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4 , 5 ])
    ax.set_xticklabels(['116_ref' , '116_mut' , '' , 'huvec_ref' , 'huvec_mut' ] , fontsize = 10)
    ax.set_ylabel('log2FC' , fontsize = 20)
    ax.set_xlabel('Type')
    ax.set_xlim((0.5 , 5.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig



hct116_data = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HCT116_deseq2_norm.csv' , header = 0)
huvec_data = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HUVEC_deseq2_norm.csv' , header = 0)



mask_ref = (hct116_data['Gene_Name'].str.replace('seq', '').astype(int) % 2) == 1
hct116_ref = hct116_data[mask_ref]
mask_mut = (hct116_data['Gene_Name'].str.replace('seq', '').astype(int) % 2) == 0
hct116_mut = hct116_data[mask_mut]


mask_ref = (huvec_data['Gene_Name'].str.replace('seq', '').astype(int) % 2) == 1
huvec_ref = huvec_data[mask_ref]
mask_mut = (huvec_data['Gene_Name'].str.replace('seq', '').astype(int) % 2) == 0
huvec_mut = huvec_data[mask_mut]




Box_plot_4cellline([hct116_ref['log2FoldChange'] , hct116_mut['log2FoldChange'] , huvec_ref['log2FoldChange'] , huvec_mut['log2FoldChange']] , -8 , 8)






























