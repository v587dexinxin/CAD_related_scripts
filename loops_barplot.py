# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 16:38:40 2024

@author: lenovo
"""



from __future__ import division
import numpy as np
import pandas as pd
import csv
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice
import os


    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()    
    
    
    
loops1 = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_control_peaks_one_anchors_binding_loops.bedpe' , header = None)
loops2 = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_ls_peaks_one_anchors_binding_loops.bedpe' , header = None)

x = [0 , 1]
y = [len(loops1) , len(loops2)]


Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6


fig = plt.figure(figsize = (12 , 12))
ax = fig.add_axes([Left  , HB , width , HH])
ax.bar(x, y , color = 'c')

ax.set_xlim((-0.5 , 1.5))
ax.set_ylim((0 , 18000))
ax.set_title('control VS LS loops number')
ax.set_xticks(x)
ax.set_xticklabels(['control' , 'LS'])
# ax.set_xlabel('')
for i in range(len(x)):
    ax.text(x[i]  , y[i] + 20 , str(y[i]))


run_Plot(fig , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\plots\\HUVEC_control_Vs_LS_Loops_barplot.pdf')




