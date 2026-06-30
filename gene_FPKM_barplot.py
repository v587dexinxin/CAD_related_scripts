# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 05:40:09 2024

@author: lenovo
"""


import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
import matplotlib.pyplot as plt


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    


union = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0 , index_col=0)
union['ST_FPKM'] = (union['ST_R1_FPKM'] + union['ST_R2_FPKM']) / 2
union['LS_FPKM'] = (union['LS_R1_FPKM'] + union['LS_R2_FPKM']) / 2
union['OS_FPKM'] = (union['OS_R1_FPKM'] + union['OS_R2_FPKM']) / 2



data = []

for n in ['KLF2' , 'KLF4' , 'NOS3' , 'BACH1' , 'BMP4']:
    genes = union.loc[n]
    data.append((n , np.log10(genes['ST_FPKM'] + 1), 'ST'))
    data.append((n , np.log10(genes['LS_FPKM'] + 1), 'LS'))
    data.append((n , np.log10(genes['OS_FPKM'] + 1), 'OS'))
    
data = pd.DataFrame(data , columns=['gene_name' , 'FPKM' , 'HUVEC_classify'])    


fig, axs = plt.subplots(2, 1,figsize=(20,4))
axs = axs.flatten() 
### contact distance
sns.barplot(x="gene_name", y="FPKM", hue='HUVEC_classify' , data=data, ax=axs[0])



# ### contact strength logq
# sns.barplot(x="PP_type", y="count", data=plotdata, ax=axs[2])

run_Plot(fig, 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\selected_genes_FPKM_barplot.pdf')


