# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 11:25:59 2024

@author: lenovo
"""

from heapq import merge
from itertools import count, islice
# from contextlib2 import ExitStack
from matplotlib.backends.backend_pdf import PdfPages
from random import random
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# from palettable.colorbrewer.qualitative import Dark2_8
import os, sys, re, time, subprocess, multiprocessing, gc, bisect, math
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
from matplotlib_venn import venn2, venn2_circles





def Volcano_Gene_Plot(fil,title,fc = 0.5 , q = 0.01):
    """
    """
    # chrom = ['chr' + str(i) for i in range(1 , 23)] + ['chrX' , 'chrY']
    f = open(fil,'r')
    Gene = []
    for line in islice(f,1,None):
        line= line.strip().split(',')
        if line[4] == 'NA' or line[8] == 'NA':
            continue
        else:
            Gene.append((line[0],line[4],-math.log10(float(line[8]) + 10**-323)))
    f.close()
    
    Gene_type = np.dtype({'names':['Gene_name' , 'FC' , 'q'],
                          'formats':['U64' , np.float64,np.float64]})
    
    Gene = np.array(Gene,dtype = Gene_type)
    # P = np.log2(Gene['M'].sum() / Gene['P'].sum())    
    
    NC_bound = fc 
    si_bound = -fc
    
    NC_mask = (Gene['FC']> NC_bound) & (Gene['q'] > -math.log10(q))
    NC_Genes = Gene[NC_mask]
    
    si_mask = (Gene['FC'] < si_bound) & (Gene['q'] > -math.log10(q))
    si_Genes = Gene[si_mask]
    
    Non_Genes = Gene[~(NC_mask | si_mask)] 
    
    
            
    fig,ax = plt.subplots(1)
    ax.scatter(NC_Genes['FC'],NC_Genes['q'], s= 10, c = 'red')
    ax.scatter(si_Genes['FC'],si_Genes['q'], s= 10, c= 'blue')
    ax.scatter(Non_Genes['FC'],Non_Genes['q'], s= 10, c = 'gray')
    ax.plot([-15,15],[-math.log10(q),-math.log10(q)], ls = '--', c = 'black', lw = 1.0)
    ax.plot([0,0],[-2,400], ls = '--', c = 'black', lw = 1.0)
    ax.plot([NC_bound,NC_bound],[-2,400], ls = '--', c = 'red', lw = 1.0)
    ax.plot([si_bound,si_bound],[-2,400],ls = '--',c = 'blue', lw = 1.0)
    
    ax.set_xticks([-8,-4,NC_bound,0,si_bound,4,8])
    ax.set_xticklabels(['-8','-4',str(NC_bound),'0',str(si_bound),'4','8'])
    ax.set_xlabel('log2FoldChange', size = 15)
    ax.set_ylabel('-log10(q-value)',size = 15)
    ax.set_ylim(-2,400)
    ax.set_xlim(-8,8)
    ax.text(3,150,'Enhancer : %d' % len(NC_Genes))
    ax.text(-6,150,'Silencer : %d' % len(si_Genes))
    ax.set_title(title)
    
    return fig,NC_Genes,si_Genes


def plot_venn2(n1 , n2 , n3 , title , out):
    fig = plt.figure(figsize = (10, 10))
    venn2(subsets=(n1 , n2 , n3), set_labels=('HCT116', 'HUVEC'))
    for text_obj in fig.findobj(matplotlib.text.Text):
        text_obj.set_fontsize(20)
        
    plt.title(title , fontsize = 20)
    run_Plot(fig , out)
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
        




pp1 = PdfPages('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HCT116_enhancer_VS_silencer_scatter_q0.05_fc0.5.pdf')


fig , enhancer_h , silencer_h = Volcano_Gene_Plot('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HCT116_deseq2_norm.csv','HCT116_Enhancer_VS_silencer',fc = 0.3, q = 0.01)

pp1.savefig(fig)
# pp1.close() 


pp1 = PdfPages('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HUVEC_enhancer_VS_silencer_scatter_q0.05_fc0.5.pdf')
fig , enhancer_hu , silencer_hu = Volcano_Gene_Plot('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HUVEC_deseq2_norm.csv','HUVEC_Enhancer_VS_silencer',fc = 0.3, q = 0.01)
pp1.savefig(fig)
pp1.close() 






#####HCT116_VS_HUVEC_Venn2_plot
enhancer_hu = list(pd.DataFrame(enhancer_hu)['Gene_name'])
enhancer_h = list(pd.DataFrame(enhancer_h)['Gene_name'])
enhancer_common = [x for x in enhancer_hu if x in enhancer_h]

silencer_hu = list(pd.DataFrame(silencer_hu)['Gene_name'])
silencer_h = list(pd.DataFrame(silencer_h)['Gene_name'])
silencer_common = [x for x in silencer_hu if x in silencer_h]


plot_venn2(len(enhancer_h) - len(enhancer_common), len(enhancer_hu) - len(enhancer_common) , len(enhancer_common) , 'Enhancer' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\DEseq2\\HCT116_HUVEC_enhancer_venn2.pdf')
plot_venn2(len(silencer_h) - len(silencer_common), len(silencer_hu) - len(silencer_common) , len(silencer_common) , 'Silencer' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\DEseq2\\HCT116_HUVEC_silencer_venn2.pdf')



#####Load_DESeq2_data
hct116_data = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HCT116_deseq2_norm.csv' , header = 0)
huvec_data = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HUVEC_deseq2_norm.csv' , header = 0)




enhancer_116 = hct116_data[hct116_data['Gene_Name'].isin(enhancer_h)]
enhancer_huvec = huvec_data[huvec_data['Gene_Name'].isin(enhancer_hu)]

silencer_116 = hct116_data[hct116_data['Gene_Name'].isin(silencer_h)]
silencer_huvec = huvec_data[huvec_data['Gene_Name'].isin(silencer_hu)]

enhancer_116.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HCT116_enhancer_deseq2_norm.csv' , header = True , index = None)
silencer_116.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HCT116_silencer_deseq2_norm.csv' , header = True , index = None)

enhancer_huvec.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HUVEC_enhancer_deseq2_norm.csv' , header = True , index = None)
silencer_huvec.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\HUVEC_silencer_deseq2_norm.csv' , header = True , index = None)










