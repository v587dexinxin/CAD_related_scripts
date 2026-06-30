# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 12:20:51 2024

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
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
def Load_peaks(file , peaks_type):
    sz = os.path.getsize(file)
    if sz != 0:
        peaks = pd.read_table(file , header = None)
        if peaks_type == 'narrow':
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue' , 'lengtn']
        else:
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue']
    else:
        peaks = []

    return peaks



def Specific_Common_peaks(peaks1 , peaks2):
    n = 0 ; speci1 = [] ; speci2 = [] ; common = []
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    
    for g in chrom:
        tmp1 = peaks1[peaks1['chr'] == g]
        tmp2 = peaks2[peaks2['chr'] == g]
        for i in tmp1.index:
            start = tmp1.loc[i]['start']
            end = tmp1.loc[i]['end']
            mask = (tmp2['start'] <= end) & (tmp2['end'] >= start)
            overlap = tmp2[mask]
            if len(overlap) != 0:
                n += 1
                common.append((g , start , end))
            else:
                speci1.append((g , start , end))
        for i in tmp2.index:
            start = tmp2.loc[i]['start']
            end = tmp2.loc[i]['end']
            mask = (tmp1['start'] <= end) & (tmp1['end'] >= start)
            overlap = tmp1[mask]
            if len(overlap) != 0:
                pass
            else:
                speci2.append((g , start , end))
                
    print (n)
    common = pd.DataFrame(common)
    common.columns = ['chr' , 'start' , 'end']
    speci1 = pd.DataFrame(speci1)
    speci1.columns = ['chr' , 'start' , 'end']
    speci2 = pd.DataFrame(speci2)
    speci2.columns = ['chr' , 'start' , 'end']
    
    return speci1 , speci2 , common

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

def Box_plot_4cellline(data , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})
    ax.boxplot(data[2] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[3] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4 , 5 ])
    ax.set_xticklabels(['K562' , 'HCT116' , '' , 'K562' , 'HCT116' ] , fontsize = 10)
    ax.set_ylabel('FPKM' , fontsize = 20)
    ax.set_xlabel('K562_specific_peaks_VS_HCT116_specific_peaks')
    ax.set_xlim((0.5 , 5.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig



def Box_plot_2cellline(data , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 ])
    ax.set_xticklabels(['HCT116_Enhancer' , 'HCT116_Silencer'] , fontsize = 10)
    ax.set_ylabel('FPKM' , fontsize = 20)
    ax.set_xlabel('HCT116_enhancer_VS_silencer_related_genes')
    ax.set_xlim((0.5 , 2.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig







##-------------0.1FA_specific_peaks-------

enhancer_116 = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\position\\HCT116_enhancer_deseq2_pos_+-1bp.bed' , header = None)
enhancer_116.columns = ['chr' , 'start' , 'end']


silencer_116 = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\position\\HCT116_silencer_deseq2_pos_+-1bp.bed' , header = None)
silencer_116.columns = ['chr' , 'start' , 'end']


##----------------RNA---------------------

RNA = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0 , sep = ',')
RNA = RNA.drop_duplicates(subset = ['Gene_Name'] , keep = 'first')
RNA['K562_FPKM'] = (RNA['K562_WT_R1_FPKM'] + RNA['K562_WT_R2_FPKM']) / 2
RNA['HCT116_FPKM'] = (RNA['HCT116_WT_R1_FPKM'] + RNA['HCT116_WT_R2_FPKM']) / 2

# expressed_rna = RNA[(RNA['K562_FPKM'] >= 2) | (RNA['HCT116_FPKM'] >= 2)]

n = 0 ; m = 0
enhancer_genes = [] 
for i in enhancer_116.index:
    g = enhancer_116.loc[i]['chr']
    pos = enhancer_116.loc[i]['start'] + 1
    tmp_RNA = RNA[RNA['Chr'] == g]
    overlap = tmp_RNA[(tmp_RNA['Start'] <= pos) & (tmp_RNA['End'] >= pos)]
    if len(overlap != 0):
        # print (g , pos)
        n += 1
        enhancer_genes.append(overlap.iloc[0]['HCT116_FPKM'])
    else:
        my_list = list(tmp_RNA['Start']) + list(tmp_RNA['End']) 
        closest_element = min(my_list, key=lambda x: abs(x - pos))
        genes = tmp_RNA[(tmp_RNA['Start'] == closest_element) | (tmp_RNA['End'] == closest_element)]
        if len(genes) == 0:
            print (g , pos)
            m += 1
        # enhancer_genes.append(genes.iloc[0]['HCT116_FPKM'])
            
            
            
silencer_genes = [] 
for i in silencer_116.index:
    g = silencer_116.loc[i]['chr']
    pos = silencer_116.loc[i]['start'] + 1
    tmp_RNA = RNA[RNA['Chr'] == g]
    overlap = tmp_RNA[(tmp_RNA['Start'] <= pos) & (tmp_RNA['End'] >= pos)]
    if len(overlap != 0):
        # print (g , pos)
        n += 1
        silencer_genes.append(overlap.iloc[0]['HCT116_FPKM'])
    else:
        my_list = list(tmp_RNA['Start']) + list(tmp_RNA['End']) 
        closest_element = min(my_list, key=lambda x: abs(x - pos))
        genes = tmp_RNA[(tmp_RNA['Start'] == closest_element) | (tmp_RNA['End'] == closest_element)]
        if len(genes) == 0:
            print (g , pos)
            m += 1
        # silencer_genes.append(genes.iloc[0]['HCT116_FPKM'])
        
        
        
        
        
        



genes = []
for i in expressed_rna.index:
    gene_name = expressed_rna.loc[i]['Gene_Name']
    g = expressed_rna.loc[i]['Chr']
    strand = expressed_rna.loc[i]['Strand']
    start = expressed_rna.loc[i]['Start']
    end = expressed_rna.loc[i]['End']
    fpkm_562 = expressed_rna.loc[i]['K562_FPKM']
    fpkm_116 = expressed_rna.loc[i]['HCT116_FPKM']
    if strand == '+':
        start = start - 2000
    else:
        end = end + 2000
    genes.append((gene_name , g , strand , start , end , fpkm_562 , fpkm_116))
    
    
genes = pd.DataFrame(genes)
genes.columns = ['Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'fpkm_562' , 'fpkm_116']

###--------------------peaks_genes_overlap------------------------

speci562_peak_genes = [] ; speci116_peak_genes = [] 
            
    
def peaks_related_genes(peaks , genes):
    peaks_genes = []
    for g in chrom:
        print (g)
        tmp_genes = genes[genes['chr'] == g]
        tmp_peaks = peaks[peaks['chr'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            mask = (tmp_genes['start'] <= end) & (tmp_genes['end'] >= start)
            overlap = tmp_genes[mask]
            if len(overlap) != 0:
                for j in overlap.index:
                    peaks_genes.append((overlap.loc[j]['Gene_name'] , overlap.loc[j]['fpkm_562'] , overlap.loc[j]['fpkm_116']))
            else:
                pass
    peaks_genes = pd.DataFrame(peaks_genes)
    peaks_genes.columns = ['Gene_name' , 'fpkm_562' , 'fpkm_116']
    peaks_genes = peaks_genes.drop_duplicates(subset = 'Gene_name')
    return peaks_genes
            


speci562_peak_genes = peaks_related_genes(speci_562 , genes)

speci116_peak_genes = peaks_related_genes(speci_116 , genes)






data = [list(speci562_peak_genes['fpkm_562']) , list(speci562_peak_genes['fpkm_116']) , list(speci116_peak_genes['fpkm_562']) , list(speci116_peak_genes['fpkm_116'])]























