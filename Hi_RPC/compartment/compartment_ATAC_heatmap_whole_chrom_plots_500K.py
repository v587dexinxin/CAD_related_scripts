# -*- coding: utf-8 -*-
"""
Created on Wed May 20 17:31:14 2026

@author: lenovo
"""

from __future__ import division
import numpy as np
import cooler
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
from matplotlib.backends.backend_pdf import PdfPages
import os
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

#---------------------------------------------------Functions-------------------------------------------------------------
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
    
def caxis_S_vertical(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'y', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
def caxis_S_horizontal(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    
def UpdateDI(DI):
    """
    """
    New_DI = []
    New_index = []

    for index in range(len(DI) - 1):
        if DI[index] * DI[index + 1] < 0:
            New_DI.append(DI[index])
            New_DI.append(0)
            New_index.append(index)
            New_index.append(index + 0.5)
        else:
            New_DI.append(DI[index])
            New_index.append(index)
    
    return np.array(New_index), np.array(New_DI)

def add_ax(sig , loc , cell):
    sig1 = [] ; sig2 = [] ; index1 = [] ; index2 = []
    for i in range(len(sig)):
        if sig[i] > 0:
            index1.append(i)
            sig1.append(sig[i])
        elif sig[i] < 0:
            index2.append(i)
            sig2.append(sig[i])
        else:
            pass
    index1 = np.array(index1)
    index2 = np.array(index2)
    sig1 = np.array(sig1)
    sig2 = np.array(sig2)

    
    ax = fig.add_axes(loc)
    ax.bar(index1 , sig1 , 1 , color = 'gold')
    ax.bar(index2 , sig2 , 1 , color = 'midnightblue')
    # ax.fill_between(PCA_index1 , PCA1 , where = PCA1 >= 0 , facecolor = 'gold' , edgecolor = 'none' )
    # ax.fill_between(PCA_index1 , PCA1 , where = PCA1 <= 0 , facecolor = 'midnightblue' , edgecolor = 'none' )
    ax.set_xlim((-0.5 , len(sig) - 0.5))
    ytick = [round(sig.min() * 0.6667 , 2) , 0.00 , round(sig.max() * 0.6667, 2)]
    ax.set_yticks(ytick)
    ax.set_ylim((sig.min() * 1.1 , sig.max()* 1.1))
    ax.set_ylabel(cell,fontsize=20,rotation = 'horizontal' , labelpad = 45)        
    print (len(sig))
    return ax


def chrom_length(file):
    hg38 = {}
    data = pd.read_table(file , header = None , sep = '\t')
    for i in range(len(data)):
        g = data.loc[i][0]
        length = data.loc[i][1]
        hg38[g] = length
    return hg38
    


def bedgraph_to_sig(file):
    data = pd.read_table(file , header = None , sep = '\t')
    data.columns = ['chr' , 'start' , 'end' , 'pca']
    pca = {}
    for g in chrom:
        print (g)
        g_length = hg38[g] // R + 1
        pca[g] = np.zeros(g_length)
        tmp = data[data['chr'] == g]
        for i in range(len(tmp)):
            start = tmp.iloc[i]['start'] // R 
            end = tmp.iloc[i]['end'] // R 
            p = tmp.iloc[i]['pca']
            for j in range(start , end):
                pca[g][j] += p
    return pca
            
    



def ATAC_bedgraph_to_sig(ATACData , R):
    
    ATACData.columns = ['chr' , 'start' , 'end' , 'score']
    sig = {}
    for g in chrom:
        print (g)
        g_length = hg38[g] // R + 1
        tmp = ATACData[ATACData['chr'] == g]
        tmp['start_res'] = tmp['start'] // R
        tmp = tmp.groupby(['start_res'], as_index=False).agg({'score': 'sum', 'start_res': 'first'})
        yindex = []
        for i in range(g_length):
            if i not in list(tmp['start_res']):
                yindex.append(i)
        new_rows = pd.DataFrame({'score':[0 for x in range(len(yindex))] , 'start_res':yindex})
        tmp = pd.concat([tmp , new_rows])
        tmp = tmp.sort_values(by = ['start_res'])

        sig[g] = tmp
        
    return sig






def Load_PC_Data(pc_file):
    pc = pd.read_table(pc_file , header = 0)
    pc = pc.fillna(0)
    return pc








#-----------------------------------------------Files-------------------------------------------------------------------------    



cell = ['WT' , 'LS' , 'OS']
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX' ]
R = 500000
res = '500K'

hg38 = chrom_length('/scratch/2026-05-18/bio-shenw/ref/Human/hg38/hg38.chrom.size')




PCFolder = '/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/compartment/'
Outfolder = '/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/compartment/plots'




###############ATAC_Data

wt_ATAC = pd.read_table('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/mapping/bam2/all_reps/signals/bedgraph/WT_ATAC.RPKM.rescaled_0.882_housekeeping_stable.bedgraph' , header = None , sep = '\t')
ls_ATAC = pd.read_table('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/mapping/bam2/all_reps/signals/bedgraph/LSS_ATAC.RPKM.rescaled_1.712_housekeeping_stable.bedgraph' , header = None , sep = '\t')
os_ATAC = pd.read_table('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/mapping/bam2/all_reps/signals/bedgraph/HUVEC_ATAC_OS_allreps_RPKM_10bp.bedgraph' , header = None , sep = '\t')



wt_atac = ATAC_bedgraph_to_sig(wt_ATAC , 100000)
ls_atac = ATAC_bedgraph_to_sig(ls_ATAC , 100000)
os_atac = ATAC_bedgraph_to_sig(os_ATAC , 100000)


ATAC_Data = {'WT' : wt_atac , 'LS' : ls_atac , 'OS' : os_atac}




##############PC_Data

wt_pc = Load_PC_Data('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/compartment/res_500k/HUVEC_WT_ATACphased_500kb.cis.vecs.tsv')
ls_pc = Load_PC_Data('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/compartment/res_500k/HUVEC_LS_ATACphased_500kb.cis.vecs.tsv')
os_pc = Load_PC_Data('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/compartment/res_500k/HUVEC_OS_ATACphased_500kb.cis.vecs.tsv')


PC_Data = {'WT' : wt_pc , 'LS' : ls_pc , 'OS' : os_pc}


#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.19 ; HB = 0.17 ; width = 0.5 ; HH = 0.5

for n in ['E1' , 'E2' , 'E3']:    
    for c in cell:
        HiCFolder = '/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool'
        HiCFil = 'HUVEC_' + c + '_merged.hg38.nodups.mapq_30.1000.mcool::/resolutions/500000'
        HiCSource = os.path.join(HiCFolder , HiCFil)
        HiCData = cooler.Cooler(HiCSource) 
        
        PCData = ATAC_bedgraph_to_sig(PC_Data[c][['chrom' , 'start' , 'end' , n]] , 500000)
        atacData = ATAC_Data[c]
        
    
        OutFil = c + '_' + res + '_Heatmap_balanced_compartment_selected_pca' + n + '.pdf'
        pp = PdfPages(os.path.join(Outfolder , OutFil))
        # g = chrom[0]
        for g in chrom:
            matrix = HiCData.matrix(balance=True).fetch(g)
            matrix[np.isnan(matrix)] = 0
            sig = PCData[g]['score']
            atac = atacData[g]['score']
            nonzero = matrix[np.nonzero(matrix)]
            vmax = np.percentile(nonzero, 95)
            vmin = nonzero.min()
            ## Heatmap Plotting
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left  , HB , width , HH])
            sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                           extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, vmin = vmin , origin = 'lower')
            cxlim = ax.get_xlim()
            cylim = ax.get_ylim()
            ## Ticks and Labels
            ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
            pos = [t * R for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels , fontsize=10)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels , fontsize=10 , rotation = 'horizontal')
            ax.set_xlabel(g , fontsize=20,labelpad=30)
    
                            
    #b3b3b3
            ax.set_xlim(cxlim)
            ax.set_ylim(cylim)                    
            ## Colorbar
            ax = fig.add_axes([Left + 0.6 , HB - 0.12 , 0.1 , 0.035])
            cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
            cbar.set_ticks([vmin , vmax])
    
            
            ##PCA Tracks
            
            ax1 = add_ax(np.array(sig) , [Left, HB + width , width , 0.1] , 'PC1')
            caxis_S_horizontal(ax1, 'black')
            
            ##ATAC Tracks
            color = 'green'
            tmp_sig = atac
            # length = int(len(matrix) * R/1000)
            # new = np.zeros(length)
            # for i in atac:
            #     start = i['start'] // 1000
            #     end = i['end'] // 1000
            #     score = i['score']
            #     new[start] += score
            vm = tmp_sig.max()
            ax = fig.add_axes([Left , HB + width + 0.1, width , 0.1])
            ax.fill_between(np.arange(len(tmp_sig)) , np.array(tmp_sig) , facecolor = color , edgecolor = 'none' )
            ax.set_xlabel('ATAC',fontsize=15)
            ax.set_ylim((0 , vm * 1.1))
            ax.set_xlim((0 , len(tmp_sig)))
            
            caxis_S_horizontal(ax1, 'black')
            pp.savefig(fig)
            plt.close(fig)
        pp.close()