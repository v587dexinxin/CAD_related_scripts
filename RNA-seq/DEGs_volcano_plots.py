# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 18:09:26 2026

@author: lenovo
"""

from heapq import merge
from itertools import count, islice
# from contextlib2 import ExitStack
from matplotlib.backends.backend_pdf import PdfPages
from random import random
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# from palettable.colorbrewer.qualitative import Dark2_8
import os, sys, re, time, subprocess, multiprocessing, gc, bisect, math
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd




def Volcano_Gene_Plot(fil,title, s1 , s2 , fc = 0.5 , q = 0.01):
    """
    """
    f = open(fil,'r')
    Gene = []
    for line in islice(f,1,None):
        line= line.strip().split(',')
        if line[6] == 'NA' or line[10] == 'NA' or (line[1] not in chrom):
            continue
        else:
            Gene.append((line[0],line[6],-math.log10(float(line[10]) + 10**-323)))
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
    ax.plot([0,0],[-2,200], ls = '--', c = 'black', lw = 1.0)
    ax.plot([NC_bound,NC_bound],[-2,200], ls = '--', c = 'red', lw = 1.0)
    ax.plot([si_bound,si_bound],[-2,200],ls = '--',c = 'blue', lw = 1.0)
    
    ax.set_xticks([-8,-4,NC_bound,0,si_bound,4,8])
    ax.set_xticklabels(['-8','-4',str(NC_bound),'0',str(si_bound),'4','8'])
    ax.set_xlabel('log2FoldChange', size = 15)
    ax.set_ylabel('-log10(q-value)',size = 15)
    ax.set_ylim(-2,200)
    ax.set_xlim(-15,15)
    ax.text(3,150,s1 + '_up_genes : %d' % len(NC_Genes))
    ax.text(-14,150,s2 + '_down_genes : %d' % len(si_Genes))
    ax.set_title(title)
    
    return fig


def Get_gene_start(genes):
    a = genes[genes['Strand'] == '+']
    b = genes[genes['Strand'] == '-']
    a['gene_start'] = a['Start']
    b['gene_start'] = b['End']
    tmp = pd.concat([a , b] , axis = 0)
    tmp = tmp.sort_values(by = ['Chr' , 'gene_start'])
    tmp['gene_start_end'] = tmp['gene_start'] + 1
    genes_start = tmp[['Chr' , 'gene_start' , 'gene_start_end' , 'Gene_Name']]
    return (genes_start)


def Get_diff_genes(genes , fc , q):
    '''
    '''
    
    wt_RNA = RNA[(RNA['Chr'].isin(chrom)) & (RNA['log2FoldChange'] >= fc) & (RNA['padj'] <= q)]
    hemin_RNA = RNA[(RNA['Chr'].isin(chrom)) & (RNA['log2FoldChange'] <= -fc) & (RNA['padj'] <= q)]
    stable_RNA = RNA[(RNA['Chr'].isin(chrom)) & ~(((RNA['log2FoldChange'] >= fc) & (RNA['padj'] <= q)) | ((RNA['log2FoldChange'] <= -fc) & (RNA['padj'] <= q)))]
    
    return wt_RNA , hemin_RNA , stable_RNA




chrom = ['chr' + str(i) for i in range(1 , 23)] + ['chrX']







###############LS_VS_OS#############


fc = 1 ; q = 0.05

RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\HUVEC_DEGs_LS_VS_OS.csv' , header = 0)

LS_up , LS_down , stable = Get_diff_genes(RNA , fc , q)

LS_up.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_up_genes_q0.05_fc0.5.csv' , header = True , index = None)
LS_down.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_down_genes_q0.05_fc0.5.csv' , header = True , index = None)


LS_up['Gene_Name'].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_up_genes_q0.05_fc0.5_gene_name.txt' , header = None , index = None , sep = '\t')
LS_down['Gene_Name'].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_down_genes_q0.05_fc0.5_gene_name.txt' , header = None , index = None , sep = '\t')



pp1 = PdfPages('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_VS_OS_RNA_DEGs_scatter_q0.05_fc0.5_1.pdf')


fig = Volcano_Gene_Plot('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\HUVEC_DEGs_LS_VS_OS.csv','LS_VS_OS_DEGs', 'LS' , 'LS' , fc, q)

pp1.savefig(fig)

pp1.close() 






###############HCT116_WT_VS_NP13#############

def Get_diff_genes_1(genes , fc , q):
    '''
    '''
    
    wt_RNA = RNA[(RNA['log2FoldChange'] >= fc) & (RNA['padj'] <= q)]
    hemin_RNA = RNA[(RNA['log2FoldChange'] <= -fc) & (RNA['padj'] <= q)]
    stable_RNA = RNA[~(((RNA['log2FoldChange'] >= fc) & (RNA['padj'] <= q)) | ((RNA['log2FoldChange'] <= -fc) & (RNA['padj'] <= q)))]
    
    return wt_RNA , hemin_RNA , stable_RNA



def Volcano_Gene_Plot_1(fil,title, s1 , s2 , fc = 0.5 , q = 0.01):
    """
    """
    f = open(fil,'r')
    Gene = []
    for line in islice(f,1,None):
        line= line.strip().split(',')
        if line[4] == 'NA' or line[-1] == 'NA':
            continue
        else:
            Gene.append((line[0],line[4],-math.log10(float(line[-1]) + 10**-323)))
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
    ax.plot([0,0],[-2,200], ls = '--', c = 'black', lw = 1.0)
    ax.plot([NC_bound,NC_bound],[-2,200], ls = '--', c = 'red', lw = 1.0)
    ax.plot([si_bound,si_bound],[-2,200],ls = '--',c = 'blue', lw = 1.0)
    
    ax.set_xticks([-8,-4,NC_bound,0,si_bound,4,8])
    ax.set_xticklabels(['-8','-4',str(NC_bound),'0',str(si_bound),'4','8'])
    ax.set_xlabel('log2FoldChange', size = 15)
    ax.set_ylabel('-log10(q-value)',size = 15)
    ax.set_ylim(-2,200)
    ax.set_xlim(-10,10)
    ax.text(3,90,s1 + '_up_genes : %d' % len(NC_Genes))
    ax.text(-9,90,s2 + '_down_genes : %d' % len(si_Genes))
    ax.set_title(title)
    
    return fig




fc = 1; q = 0.05

RNA = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\HCT116_WT_VS_NP13_KO_deseq2.csv' , header = 0)

LS_down , LS_up , stable = Get_diff_genes_1(RNA , fc , q)

LS_up.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\KO13_up_genes_q0.05_fc1.csv' , header = True , index = None)
LS_down.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\KO13_down_genes_q0.05_fc1.csv' , header = True , index = None)


LS_up['Gene_Name'].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\KO13_up_genes_q0.05_fc1_gene_name.txt' , header = None , index = None , sep = '\t')
LS_down['Gene_Name'].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\KO13_down_genes_q0.05_fc1_gene_name.txt' , header = None , index = None , sep = '\t')



pp1 = PdfPages('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\WT_VS_KO13_RNA_DEGs_scatter_q0.05_fc1.pdf')


fig = Volcano_Gene_Plot_1('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\HCT116_WT_VS_NP13_KO_deseq2.csv', 'WT_VS_KO13' , 'KO' , 'KO' , fc = 1, q = 0.05)

pp1.savefig(fig)

pp1.close() 









