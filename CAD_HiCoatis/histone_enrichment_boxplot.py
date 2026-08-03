# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:56:57 2026

@author: lenovo
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import pyBigWig
import seaborn as sns
from scipy.interpolate import  interp1d
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy
from scipy import stats
from scipy.stats import ttest_ind
import pandas as pd

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')


    
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['U8' , np.int64 , np.int64]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int64 , np.int64 , np.float64]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , '20' , '21' , '22' , 'X']
res = 200000



def Sig_To_1K(signal):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 1000 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 1000
            New_Data[g][start] += line['value']
    
    return New_Data

def Get_random(Chip):
    r = []
    for i in Chip:
        for j in Chip[i]:
            r.append(j)
    return(r)





def Get_signal_cluster(chip , pc_data):
    cluster = pc_data.keys()
    sig = {} 
    for c in cluster:
        sig[c] = []
        for i in pc_data[c].index:
            g = pc_data[c].loc[i]['chr'].lstrip('chr')
            start = pc_data[c].loc[i]['start'] // 1000
            end = pc_data[c].loc[i]['end'] // 1000
            if start == end:
                sig_mean = chip[g][start]
            else:
                sig_chip = chip[g][start:end]
                sig_mean = np.mean(sig_chip)
            sig[c].append(sig_mean)
    # sig['All'] = [sig['All'] += sig[x] for x in cluster]
    return sig


def Box_plot(data , histone):
        
    pp = PdfPages('/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/plots/Histone_enrichment_boxplot/Cardiovascular_disease_related_SNPs_' + histone + '_enrichment_boxplot.pdf')
           

    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    
    color = ['#E31A1C', '#1F78B4', '#6A3D9A']
    
    
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data['promoter'] , 
            positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[0],'linewidth':1},
            medianprops={'color':color[0],'linewidth':1},
            capprops={'color':color[0],'linewidth':1},
            whiskerprops={'color':color[0],'linewidth':1, 'linestyle':'--'})
    ax.boxplot(data['Exon'] , 
            positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[1],'linewidth':1},
            medianprops={'color':color[1],'linewidth':1},
            capprops={'color':color[1],'linewidth':1},
            whiskerprops={'color':color[1],'linewidth':1, 'linestyle':'--'})
    ax.boxplot(data['Noncoding'] , 
            positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[2],'linewidth':1},
            medianprops={'color':color[2],'linewidth':1},
            capprops={'color':color[2],'linewidth':1},
            whiskerprops={'color':color[2],'linewidth':1, 'linestyle':'--'})
                        
    d1 = np.round(ttest_ind(data['promoter'],data['Noncoding'])[1] , 5)
    d2 = np.round(ttest_ind(data['Exon'],data['Noncoding'])[1] , 5)
    
    d1 = scipy.stats.ranksums(data['promoter'] , data['Noncoding'])[1]
    d1 = scipy.stats.ranksums(data['Exon'] , data['Noncoding'])[1]
    
    ax.set_xticks([0,1,2,3,4])
    ax.set_xticklabels(['','Promoter' , 'Exon' , 'Noncoding',''] ,fontsize = 20)
    ax.set_xlabel(['d1:' + str(d1) + ',d2:' + str(d2)])
    # ax.set_ylim(-5 , 5)
    ax.set_ylabel(histone + ' signal intensity around SNPs')
    
    pp.savefig(fig)
    pp.close()             
            
            


def Box_plot_1(sig1, sig2 , histone):
        
    pp = PdfPages('/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/plots/Histone_enrichment_boxplot/HUVEC_' + histone + '_enrichment_boxplot.pdf')
           

    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    
    color = ["#BDBDBD" , '#1F78B4', "#54A24B" , "#D65F9E"]
    
    
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(sig1 , 
            positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[0],'linewidth':2},
            medianprops={'color':color[0],'linewidth':2},
            capprops={'color':color[0],'linewidth':2},
            whiskerprops={'color':color[0],'linewidth':2, 'linestyle':'--'})
    ax.boxplot(sig2 , 
            positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[1],'linewidth':2},
            medianprops={'color':color[1],'linewidth':2},
            capprops={'color':color[1],'linewidth':2},
            whiskerprops={'color':color[1],'linewidth':2, 'linestyle':'--'})
                        
    # d1 = np.round(ttest_ind(data['promoter'],data['Noncoding'])[1] , 5)

    
    d1 = scipy.stats.ranksums(sig1 , sig2)[1]

    
    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels(['','Random' , 'CAD_SNPS_region',''] ,fontsize = 20)
    ax.set_xlabel(['d1:' + str(d1) ])
    # ax.set_ylim(-5 , 5)
    ax.set_ylabel(histone + ' signal intensity around SNPs' ,fontsize = 20)
    
    pp.savefig(fig)
    pp.close()             

    



def Box_plot_2(data , histone):
        
    pp = PdfPages('/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/plots/Histone_enrichment_boxplot/HUVEC_' + histone + '_enrichment_boxplot.pdf')
           

    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    
    color = ["#BDBDBD" , '#1F78B4', "#54A24B" , "#D65F9E"]
    
    
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , 
            positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[0],'linewidth':2},
            medianprops={'color':color[0],'linewidth':2},
            capprops={'color':color[0],'linewidth':2},
            whiskerprops={'color':color[0],'linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[1] , 
            positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[1],'linewidth':2},
            medianprops={'color':color[1],'linewidth':2},
            capprops={'color':color[1],'linewidth':2},
            whiskerprops={'color':color[1],'linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[2] , 
            positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[2],'linewidth':2},
            medianprops={'color':color[2],'linewidth':2},
            capprops={'color':color[2],'linewidth':2},
            whiskerprops={'color':color[2],'linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[3] , 
            positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': color[3],'linewidth':2},
            medianprops={'color':color[3],'linewidth':2},
            capprops={'color':color[3],'linewidth':2},
            whiskerprops={'color':color[3],'linewidth':2, 'linestyle':'--'})
                        
    # d1 = np.round(ttest_ind(data['promoter'],data['Noncoding'])[1] , 5)

    
    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
    d3 = scipy.stats.ranksums(data[0] , data[3])[1]
    d4 = scipy.stats.ranksums(data[1] , data[2])[1]
    d5 = scipy.stats.ranksums(data[1] , data[3])[1]
    d6 = scipy.stats.ranksums(data[2] , data[3])[1]
    
    

    
    ax.set_xticks([0,1,2,3,4,5])
    ax.set_xticklabels(['','Random' , 'WT' , 'LSS' , 'OSS' , ''] ,fontsize = 20)
    ax.set_xlabel(['12:' + str(d1) + ';' + '13:' + str(d2) + ';' + '14:' + str(d3) + ';' + '23:' + str(d4) + ';' + '24:' + str(d5) + ';' + '34:' + str(d6)])
    # ax.set_ylim(-5 , 5)
    ax.set_ylabel(histone + ' signal intensity around SNPs' ,fontsize = 20)
    
    pp.savefig(fig)
    pp.close()             







    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    





                    
pc1 = pd.read_table('/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/CAD_SNPs/Hapmap/CAD_related_SNPs_LD0.8_all_+-200bp.bed' , header = None)
pc1.columns = ['chr' , 'start' , 'end']
pc_data = {'All' : pc1}
# chip1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3_new/mapping/signals/CCs_H3K4me3_R2.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3/mapping/Input.bw")


chip1 = pyBigWig.open("/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/bam2/all_reps/signals/WT.RPKM.rescaled_0.605_stable_housekeeping.bw")
chip2 = pyBigWig.open("/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/ATAC-seq/mapping/bam2/all_reps/signals/WT_ATAC.RPKM.rescaled_0.882_housekeeping_stable.bw")
chip3 = pyBigWig.open("/scratch/2026-05-25/bio-shenw/literature/HUVEC/Human_endothelial_cell_of_umbilical_vein_H3K27ac_hg38_ENCFF955PAU.bigWig")
chip4 = pyBigWig.open("/scratch/2026-05-25/bio-shenw/literature/HUVEC/Human_endothelial_cell_of_umbilical_vein_H3K4me3_hg38_ENCFF161GMO.bigWig")

chip5 = pyBigWig.open("/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/bam2/all_reps/signals/LSS.RPKM.rescaled_1.518_stable_housekeeping.bw")
chip6 = pyBigWig.open("/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/bam2/all_reps/signals/HiRPC_OS_allreps_RPKM_10bp.bw")



Chip1 = Sig_To_1K(chip1)
Chip2 = Sig_To_1K(chip2)
Chip3 = Sig_To_1K(chip3)
Chip4 = Sig_To_1K(chip4)

Chip5 = Sig_To_1K(chip5)
Chip6 = Sig_To_1K(chip6)


r1 =  Get_random(Chip1)
r2 =  Get_random(Chip2)
r3 =  Get_random(Chip3)
r4 =  Get_random(Chip4)




sig1 = Get_signal_cluster(Chip1 , pc_data)
sig2 = Get_signal_cluster(Chip2 , pc_data)
sig3 = Get_signal_cluster(Chip3 , pc_data)
sig4 = Get_signal_cluster(Chip4 , pc_data)

sig5 = Get_signal_cluster(Chip5 , pc_data)
sig6 = Get_signal_cluster(Chip6 , pc_data)






Box_plot_1(r1, sig1['All'] , 'Hi-Coatis')
Box_plot_1(r2, sig2['All'] , 'ATAC')
Box_plot_1(r3, sig3['All'] , 'H3K27ac')
Box_plot_1(r4, sig4['All'] , 'H3K4me3')



Box_plot_2([r1, sig1['All'] , sig5['All'] , sig6['All']] , 'WT_LS_OS_Hi-Coatis')

