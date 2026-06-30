# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:02:44 2023

@author: 86182
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





def Get_signal_cluster(chip):
    cluster = ['promoter' , 'Exon' , 'Noncoding']
    sig = {'promoter' : [] , 'Exon' : [] , 'Noncoding' : [] }
    for c in cluster:
        for i in pc_data[c]:
            g = i['chr'].lstrip('chr')
            start = i['start'] // 1000
            end = i['end'] // 1000
            if start == end:
                sig_mean = chip[g][start]
            else:
                sig_chip = chip[g][start:end]
                sig_mean = np.mean(sig_chip)
            sig[c].append(sig_mean)
    sig['All'] = sig['promoter'] + sig['Exon'] + sig['Noncoding']
    return sig


def Box_plot(data , histone):
        
    pp = PdfPages('/scratch/2023-01-09/bio-shenw/Lixinxin/Cardiovascular_disease/plot/Cardiovascular_disease_related_SNPs_' + histone + '_enrichment_boxplot.pdf')
           

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
        
    pp = PdfPages('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/plots/Cardiovascular_disease_related_SNPs_' + histone + '_enrichment_boxplot_2.pdf')
           

    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    
    color = ['#E31A1C', '#1F78B4', '#6A3D9A']
    
    
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
    ax.set_xticklabels(['','Random' , 'All_CAD_SNPs',''] ,fontsize = 20)
    ax.set_xlabel(['d1:' + str(d1) ])
    # ax.set_ylim(-5 , 5)
    ax.set_ylabel(histone + ' signal intensity around SNPs' ,fontsize = 20)
    
    pp.savefig(fig)
    pp.close()             

    
    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    





                    
pc1 = np.loadtxt('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/Cardiovascular_related_SNPs_Promoter_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type , encoding = 'UTF-8')
pc2 = np.loadtxt('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/Cardiovascular_related_SNPs_Exon_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/Cardiovascular_related_SNPs_Noncoding_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )

pc_data = {'promoter' : pc1 , 'Exon' : pc2 , 'Noncoding' : pc3}

# chip1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3_new/mapping/signals/CCs_H3K4me3_R2.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3/mapping/Input.bw")


chip1 = pyBigWig.open("/scratch/2024-03-11/bio-shenw/literature/Heart/Heart_left_ventricle-H3K27ac-hg38-ENCFF969BAO.bigWig")
chip2 = pyBigWig.open("/scratch/2024-03-11/bio-shenw/literature/Heart/Heart-H3K4me1-hg38-ENCFF366DMT.bigWig")
chip3 = pyBigWig.open("/scratch/2024-03-11/bio-shenw/literature/Heart/Heart-H3K4me3-hg38-ENCFF513OZY.bigWig")

Chip1 = Sig_To_1K(chip1)
Chip2 = Sig_To_1K(chip2)
Chip3 = Sig_To_1K(chip3)

r1 =  Get_random(Chip1)
r2 =  Get_random(Chip2)
r3 =  Get_random(Chip3)





sig1 = Get_signal_cluster(Chip1)
sig2 = Get_signal_cluster(Chip2)
sig3 = Get_signal_cluster(Chip3)

Box_plot_1(r1, sig1['All'] , 'H3K27ac')
Box_plot_1(r2, sig2['All'] , 'H3K4me1')
Box_plot_1(r3, sig3['All'] , 'H3K4me3')




####################################################################

        
pc1 = np.loadtxt('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/Cardiovascular_related_SNPs_Promoter_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type , encoding = 'UTF-8')
pc2 = np.loadtxt('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/Cardiovascular_related_SNPs_Exon_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/scratch/2024-03-11/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/Cardiovascular_related_SNPs_Noncoding_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )

pc_data = {'promoter' : pc1 , 'Exon' : pc2 , 'Noncoding' : pc3}

# chip1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3_new/mapping/signals/CCs_H3K4me3_R2.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3/mapping/Input.bw")


chip1 = pyBigWig.open("/scratch/2024-03-11/bio-shenw/literature/Heart/Heart-H3K36me3-hg38-ENCFF324WCU.bigwig")
chip2 = pyBigWig.open("/scratch/2024-03-11/bio-shenw/literature/Heart/Heart-H3K27me3-hg38-ENCFF055FJT.bigWig")
chip3 = pyBigWig.open("/scratch/2024-03-11/bio-shenw/literature/Heart/Heart-H3K9me3-hg38-ENCFF977PJI.bigWig")

Chip1 = Sig_To_1K(chip1)
Chip2 = Sig_To_1K(chip2)
Chip3 = Sig_To_1K(chip3)

r1 =  Get_random(Chip1)
r2 =  Get_random(Chip2)
r3 =  Get_random(Chip3)





sig1 = Get_signal_cluster(Chip1)
sig2 = Get_signal_cluster(Chip2)
sig3 = Get_signal_cluster(Chip3)

Box_plot_1(r1, sig1['All'] , 'H3K36me3')
Box_plot_1(r2, sig2['All'] , 'H3K27me3')
Box_plot_1(r3, sig3['All'] , 'H3K9me3')







