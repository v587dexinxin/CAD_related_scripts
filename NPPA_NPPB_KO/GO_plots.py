# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 22:22:16 2025

@author: lenovo
"""

import xlrd
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_datas(goFil , lanes = []):
    goData = pd.read_csv(goFil , header = 0)
    goData['-log10 (PValue)'] = -np.log10(goData['pvalue'])
    goData = goData[['Description' , '-log10 (PValue)']]
    goData = goData.loc[lanes]
    goData = goData.sort_values(by = "-log10 (PValue)")
    x = [i for i in goData['Description']]
    y = [i for i in goData['-log10 (PValue)']]
    
    return x , y
    
def Bar_plot(x , y , color):
    left, bottom, width, height = 0.45, 0.1, 0.5, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (30, 12))
    ax = fig.add_axes(size_axes)
    ax.barh(range(len(x)) , y , color = color)
    
    ax.set_yticks(range(len(x)))
    ax.set_yticklabels(x,fontsize = 30)
    ax.set_xlabel('-log10(p value)' , fontsize = 20 )

    return fig
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

    
#GO_Files
goFil_WT = 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\WT_VS_NP13_down_GO_enrichment_results.csv'
goFil_H = 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\WT_VS_NP13_up_GO_enrichment_results.csv'

x_wt , y_wt = get_datas(goFil_WT , lanes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,22])
x_H , y_H = get_datas(goFil_H , lanes = [0,4,5,6,7,8,9,10,11,13,15,33,35,41,42])



fig1 = Bar_plot(x_wt , y_wt , 'deepskyblue')
fig2 = Bar_plot(x_H , y_H , 'palevioletred')


run_Plot(fig1 , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\HCT116_WT_specific_gene_GO_barplot.pdf')
run_Plot(fig2 , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\HCT116_NP_KO13_specific_gene_GO_barplot.pdf')







