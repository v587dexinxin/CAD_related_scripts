# -*- coding: utf-8 -*-
"""
Created on Sat May 23 17:39:19 2026

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

speci_ls = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\overlaped_with_ATAC\\HUVEC_LS_specific_peaks_DiffBind_DESeq2_q0.05_fc0.5_overlaped_ATAC_Anno.csv' , header = 0)
speci_os = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\overlaped_with_ATAC\\HUVEC_OS_specific_peaks_DiffBind_DESeq2_q0.05_fc0.5_overlaped_ATAC_Anno.csv' , header = 0)

speci_ls_pro = speci_ls[speci_ls["annotation"].str.contains("Promoter", na=False)]
speci_os_pro = speci_os[speci_os["annotation"].str.contains("Promoter", na=False)]


# speci_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\bam\\DESeq2_only\\HiCoatis_LSS_vs_OSS_DESeq2_OSS_speci.bed' , header = None)
# speci_os.columns = ['seqnames' , 'start' , 'end']


deseq2_pos = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header=0)

deseq2_pos.columns = ['peak_id', 'Chr', 'Start', 'End', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'direction']

diff = deseq2_pos[deseq2_pos['padj'] <= 0.01]

oss_speci_1 = diff[diff['log2FoldChange'] >= 1]


speci_os = oss_speci_1[['Chr' , 'Start' , 'End']]
speci_os.columns = ['seqnames' , 'start' , 'end']

##----------------RNA---------------------

RNA = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0 , sep = ',')
RNA = RNA.drop_duplicates(subset = ['Gene_Name'] , keep = 'first')
RNA = RNA[~RNA["Gene_Name"].str.contains("STRG", na=False)]
RNA['LS_FPKM'] = (RNA['LS_R1_FPKM'] + RNA['LS_R2_FPKM']) / 2
RNA['OS_FPKM'] = (RNA['OS_R1_FPKM'] + RNA['OS_R2_FPKM']) / 2

expressed_rna = RNA[(RNA['LS_FPKM'] >= 2) | (RNA['OS_FPKM'] >= 2)]


ls_speci_genes = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_up_genes_q0.05_fc0.5.csv' , header = 0)
os_speci_genes = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\DEGs\\LS_down_genes_q0.05_fc0.5.csv' , header = 0)

DEGs = pd.concat([ls_speci_genes , os_speci_genes])

DEGs['LS_FPKM'] = (DEGs['LS_R1_FPKM'] + DEGs['LS_R2_FPKM']) / 2
DEGs['OS_FPKM'] = (DEGs['OS_R1_FPKM'] + DEGs['OS_R2_FPKM']) / 2


DEGs = DEGs.reset_index(drop=True)

expressed_rna = DEGs[(DEGs['LS_FPKM'] >= 2) | (DEGs['OS_FPKM'] >= 2)]

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
        Start = start - 2000
        End = start + 500
    else:
        Start = end - 500
        End = end + 2000
    genes.append((gene_name , g , strand , Start , End , fpkm_562 , fpkm_116))
    
    
genes = pd.DataFrame(genes)
genes.columns = ['Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'fpkm_ls' , 'fpkm_os']

###--------------------peaks_genes_overlap------------------------
            
    
def peaks_related_genes(peaks , genes, speci):
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
            if len(overlap) == 1:
                peaks_genes.append((overlap.iloc[0]['Gene_name'] , overlap.iloc[0]['fpkm_ls'] , overlap.iloc[0]['fpkm_os']))
            elif len(overlap) > 1:
                if speci == 'lss_speci':
                    overlap['diff'] = overlap['fpkm_ls'] / overlap['fpkm_os']
                    index_max = overlap["diff"].idxmax()
                    peaks_genes.append((overlap.loc[index_max]['Gene_name'] , overlap.loc[index_max]['fpkm_ls'] , overlap.loc[index_max]['fpkm_os']))
                elif speci == 'oss_speci':
                    overlap['diff'] = overlap['fpkm_os'] / overlap['fpkm_ls']
                    index_max = overlap["diff"].idxmax()
                    peaks_genes.append((overlap.loc[index_max]['Gene_name'] , overlap.loc[index_max]['fpkm_ls'] , overlap.loc[index_max]['fpkm_os']))
                    if len(overlap) > 1:
                        print (overlap , (overlap.loc[index_max]['Gene_name'] , overlap.loc[index_max]['fpkm_ls'] , overlap.loc[index_max]['fpkm_os']))
                    
            else:
                pass
    peaks_genes = pd.DataFrame(peaks_genes)
    peaks_genes.columns = ['Gene_name' , 'fpkm_ls' , 'fpkm_os']
    peaks_genes = peaks_genes.drop_duplicates(subset = 'Gene_name')
    return peaks_genes
            


# specils_peak_genes = peaks_related_genes(speci_ls_pro , genes)

# specios_peak_genes = peaks_related_genes(speci_os_pro , genes)


specils_peak_genes = peaks_related_genes(speci_ls , genes , 'lss_speci')

specios_peak_genes = peaks_related_genes(speci_os , genes , 'oss_speci')




data = [list(specils_peak_genes['fpkm_ls']) , list(specils_peak_genes['fpkm_os']) , list(specios_peak_genes['fpkm_ls']) , list(specios_peak_genes['fpkm_os'])]

fig = Box_plot_4cellline(data , 0 , 80)

run_Plot(fig , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig3\\Diff_peaks_related_DEGs_FPKM.pdf')

