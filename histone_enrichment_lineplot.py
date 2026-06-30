# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 22:24:45 2023

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
from scipy.interpolate import make_interp_spline

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



def Sig_To_100bp(signal):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 100 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 100
            New_Data[g][start] += line['value']
    
    return New_Data




def Get_random(Chip):
    r = []
    for i in Chip:
        for j in Chip[i]:
            r.append(j)
    return(r)





def Get_signals(chip):
    sig = [0 for x in range(60)]
    for i in pc:
        g = i['chr'].lstrip('chr')
        pos = i['start'] // 100
        start = pos - 30
        end = pos + 30
        sig_chip = chip[g][start:end].tolist() + [0] * 60
        sig_chip = sig_chip[:60]
        sig = [sig[x] + sig_chip[x] for x in range(len(sig))]
    sig = [sig[x] / len(pc) for x in range(len(sig))]
        
    return sig



            
            


def Line_plot(sig1, sig2 , histone):
        
    pp = PdfPages('/scratch/2023-01-17/bio-shenw/Lixinxin/Cardiovascular_disease/plot/Cardiovascular_disease_related_SNPs_' + histone + '_enrichment_lineplot.pdf')
           

    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    
    color = ['#E31A1C', '#1F78B4', '#6A3D9A']
    x = np.arange(0,len(sig1))
    y1 = np.array(sig1)
    y2 = np.array(sig2)
    ymax = max([max(y1) , max(y2)]) * 1.5 
    x_smooth = np.linspace(x.min() , x.max() , 400)
    y1_smooth = make_interp_spline(x, y1)(x_smooth)
    y2_smooth = make_interp_spline(x, y2)(x_smooth)
    
    
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.plot(x_smooth , y1_smooth , c = 'b' , linewidth = 1)
    ax.plot(x_smooth , y2_smooth , c = color[0] , linewidth = 1)


    
    ax.set_xticks([0,30,60])
    ax.set_xticklabels(['-3Kb','SNP','3Kb'] ,fontsize = 20)
    ax.set_xlabel('Diatance around CVD related SNPs' , fontsize = 20)
    ax.set_xlim(-1 , 61)
    ax.set_ylim(0 , ymax)
    ax.set_ylabel(histone + ' signal intensity ' ,fontsize = 20)
    
    pp.savefig(fig)
    pp.close()             

    
    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
def splitArrayByStep(lst, step):
    return [lst[i:i+step] for i in range(0, len(lst), step)]

def BGSigValue(Sigs):
    """
    计算对照组的信号值
    """
    splitSigs = splitArrayByStep(Sigs, 60)
    tempSigs2d = np.array([i for i in splitSigs if len(i) == 60])
    averageSigs2d = np.mean(tempSigs2d, axis = 0)
    return averageSigs2d
    
    
    


                    
pc = np.loadtxt('/scratch/2023-01-17/bio-shenw/Lixinxin/Cardiovascular_disease/Cardiovascular_related_SNPs_Noncoding_1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )



# chip1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3_new/mapping/signals/CCs_H3K4me3_R2.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3/mapping/Input.bw")


chip1 = pyBigWig.open("/scratch/2023-01-17/bio-shenw/Lixinxin/literature/Encode/Heart_left_ventricle-H3K27ac-hg38-ENCFF969BAO.bigWig")
chip2 = pyBigWig.open("/scratch/2023-01-17/bio-shenw/Lixinxin/literature/Encode/Heart-H3K4me1-hg38-ENCFF366DMT.bigWig")
chip3 = pyBigWig.open("/scratch/2023-01-17/bio-shenw/Lixinxin/literature/Encode/Heart-H3K4me3-hg38-ENCFF513OZY.bigWig")

Chip1 = Sig_To_100bp(chip1)
Chip2 = Sig_To_100bp(chip2)
Chip3 = Sig_To_100bp(chip3)

r1 =  Get_random(Chip1)
r2 =  Get_random(Chip2)
r3 =  Get_random(Chip3)

random1 = BGSigValue(r1)
random2 = BGSigValue(r2)
random3 = BGSigValue(r3)



sig1 = Get_signals(Chip1)
sig2 = Get_signals(Chip2)
sig3 = Get_signals(Chip3)

Line_plot(random1 , sig1 , 'H3K27ac')
Line_plot(random2 , sig2 , 'H3K4me1')
Line_plot(random3 , sig3 , 'H3K4me3')






