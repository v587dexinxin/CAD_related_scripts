# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 19:58:10 2022

@author: 86182
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

data = pd.read_table('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_related_SNPs_Promoter.txt' , header = 0)



d = []
for i in data.index:
    if data.iloc[i]['geneStrand'] == 1:
        d.append(data.iloc[i]['geneStart'] - data.iloc[i]['start'])
    elif data.iloc[i]['geneStrand'] == 2:
        d.append(data.iloc[i]['end'] - data.iloc[i]['geneEnd'])
    else:
        print (i)
        
        

a = plt.hist(d , bins=20)
x = []
for i in range(len(a[1]) - 1):
    dis = (a[1][i + 1] - a[1][i]) / 2
    x.append(dis + a[1][i])
    
y = a[0]

left, bottom, width, height = 0.35, 0.1, 0.5, 0.8
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12 , 10))
ax = fig.add_axes(size_axes)
plt.plot(x , y)    
plt.xlim(-500, 2500)
# plt.ylim(0, 300)
plt.xlabel('Distance', fontsize = 20)
plt.ylabel('SNP numbers', fontsize = 20)
plt.xticks([0 , 1000 , 2000])
ax.set_xticklabels(['TSS' , '1Kb' , '2Kb'])

plt.title('Fragment length frequency distribution', fontsize = 25)
   
run_Plot(fig , 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\Fragment length frequency distribution.pdf')    














