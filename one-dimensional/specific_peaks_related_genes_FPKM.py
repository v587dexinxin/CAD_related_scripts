# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:32:22 2026

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



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
def Load_diff_peaks(file):
    diff_peaks = pd.read_csv(file , header = 0 )
    diff_peaks = diff_peaks[diff_peaks['FDR'] <= 0.05]
    speci_562 = speci_562 = diff_peaks[diff_peaks['Fold'] < 0]
    speci_116 = diff_peaks[diff_peaks['Fold'] > 0]

    return speci_562 , speci_116






def Box_plot_4cellline(data , vmin , vmax):                
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
    ax.boxplot(data[2] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[3] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[2] , data[3])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    d2 = np.round(scipy.stats.ranksums(data[2] , data[3])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4 , 5 ])
    ax.set_xticklabels(['LS' , 'OS' , '' , 'LS' , 'OS'] , fontsize = 10)
    ax.set_ylabel('FPKM' , fontsize = 20)
    ax.set_xlabel('LS_specific_peaks:' + str(d1) + '_VS_OS_specific_peaks:' + str(d2))
    ax.set_xlim((0.5 , 5.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig


                


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

##-------------0.1FA_diff_peaks-------

speci_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\overlaped_with_ATAC\\LSS_DiffBind_DEseq2_specific_peaks_overlaped_ATAC_peaks.bed' , header = None)
speci_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\overlaped_with_ATAC\\LSS_DiffBind_DEseq2_specific_peaks_overlaped_ATAC_peaks.bed' , header = None)

speci_ls.columns = ['chr' , 'start' , 'end']
speci_os.columns = ['chr' , 'start' , 'end']

##----------------RNA---------------------

RNA = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0 , sep = ',')
RNA = RNA.drop_duplicates(subset = ['Gene_Name'] , keep = 'first')
RNA['LS_FPKM'] = (RNA['LS_R1_FPKM'] + RNA['LS_R2_FPKM']) / 2
RNA['OS_FPKM'] = (RNA['OS_R1_FPKM'] + RNA['OS_R2_FPKM']) / 2

expressed_rna = RNA[(RNA['LS_FPKM'] >= 2) | (RNA['OS_FPKM'] >= 2)]



genes = []
for i in expressed_rna.index:
    gene_name = expressed_rna.loc[i]['Gene_Name']
    g = expressed_rna.loc[i]['Chr']
    strand = expressed_rna.loc[i]['Strand']
    start = expressed_rna.loc[i]['Start']
    end = expressed_rna.loc[i]['End']
    fpkm_562 = expressed_rna.loc[i]['LS_FPKM']
    fpkm_116 = expressed_rna.loc[i]['OS_FPKM']
    if strand == '+':
        start = start - 2000
        end = start + 500
    else:
        start = end - 500
        end = end + 2000
    genes.append((gene_name , g , strand , start , end , fpkm_562 , fpkm_116))
    
    
genes = pd.DataFrame(genes)
genes.columns = ['Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'fpkm_ls' , 'fpkm_os']

###--------------------peaks_genes_overlap------------------------
            
    
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
                    peaks_genes.append((overlap.loc[j]['Gene_name'] , overlap.loc[j]['fpkm_ls'] , overlap.loc[j]['fpkm_os']))
            else:
                pass
    peaks_genes = pd.DataFrame(peaks_genes)
    peaks_genes.columns = ['Gene_name' , 'fpkm_ls' , 'fpkm_os']
    peaks_genes = peaks_genes.drop_duplicates(subset = 'Gene_name')
    return peaks_genes
            


specils_peak_genes = peaks_related_genes(speci_ls , genes)

specios_peak_genes = peaks_related_genes(speci_os , genes)






data = [list(specils_peak_genes['fpkm_ls']) , list(specils_peak_genes['fpkm_os']) , list(specios_peak_genes['fpkm_ls']) , list(specios_peak_genes['fpkm_os'])]

fig = Box_plot_4cellline(data , 0 , 20)

run_Plot(fig , 'H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_HiRPC_0.1FA\\DiffBind\\Diff_peaks_related_genes_FPKM_new.pdf')




###-------------DEGs--------

DEGs = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\DEGs\\DEGs_K562_VS_HCT116.csv' , header = 0 )
DEGs = DEGs[(DEGs['padj'] <= 0.01) & ((DEGs['log2FoldChange'] >= 1) | (DEGs['log2FoldChange'] <= -1))]

speci_562_genes = DEGs[DEGs['log2FoldChange'] > 0]

speci_116_genes = DEGs[DEGs['log2FoldChange'] < 0]


######peaks_related_genes_all_peaks

def peaks_related_genes(peaks , genes):
    peaks_genes = []
    for g in chrom:
        print (g)
        tmp_genes = genes[genes['chr'] == g]
        tmp_peaks = peaks[peaks['seqnames'] == g]
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
    # peaks_genes = peaks_genes.drop_duplicates(subset = 'Gene_name')
    return peaks_genes
            


speci562_peak_genes = peaks_related_genes(speci_562 , genes)

speci116_peak_genes = peaks_related_genes(speci_116 , genes)


common_562 = list(set(speci562_peak_genes['Gene_name']) & set(speci_562_genes['Gene_Name']))
common_116 = list(set(speci116_peak_genes['Gene_name']) & set(speci_116_genes['Gene_Name']))



precent = []
n1 = 0
for i in speci562_peak_genes.index:
    if speci562_peak_genes.loc[i]['Gene_name'] in list(speci_562_genes['Gene_Name']):
        n1 += 1
        
precent.append(n1 / len(speci562_peak_genes))


n2 = 0
for i in speci116_peak_genes.index:
    if speci116_peak_genes.loc[i]['Gene_name'] in list(speci_116_genes['Gene_Name']):
        n2 += 1
        
precent.append(n2 / len(speci116_peak_genes))



len(speci_562_genes['Gene_Name']) - len(common_562)


    









