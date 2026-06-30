# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 20:34:23 2022

@author: 86182
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



chro = [str(x) for x in range(1, 23)] + ['X']

SNPs = {}
for g in chro:
    print (g)
    out = 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed'
    data_1 = pd.read_table('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\ld_chr' + g + '_CEU.txt' , header = None , sep = ' ')
    data_1.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    data_2 = pd.read_table('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\SNP_bed\\bed_chr_' + g + '.bed' , header = None , skiprows=1 , index_col=3)
    data_2.columns = ['chr' , 'start' , 'end' , 'num' , 'strand']
    rs = list(data_1.rs1) + list(data_1.rs2)
    rs = set(rs)
    index = data_2.index
    rs_new = []
    for i in rs:
        if i not in index:
            print (i)
        else:
            rs_new.append(i)
    data_1 = data_1[data_1.rs1.isin(rs_new)]
    data_1 = data_1[data_1.rs2.isin(rs_new)]
    data_2 = data_2.loc[rs_new]
    data1_dict = {}
    for i in data_2.index:
        data1_dict[i] = data_2.loc[i].end    
    data_1.pos1 = data_1.rs1.apply(lambda x:data1_dict[x])
    data_1.pos2 = data_1.rs2.apply(lambda x:data1_dict[x])
    data_1.to_csv(out , sep='\t', index=False)
    



