# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 16:16:26 2022

@author: 86182
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    

'''
Read EBI table
'''

data = pd.read_table('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv' , header = 0 , dtype={"REPLICATION SAMPLE SIZE" : str, "CHR_ID" : str , "CHR_POS" : str , "SNP_ID_CURRENT" : str})


'''
Select disease-related data
'''
number = []

for i in range(len(data)):
    for j in data.loc[i]:
        j = str(j)
        a = j.lower()
        if ('heart' in a) or ('cardiovascular' in a):
            # print (i)
            number.append(i)
            break



data_new = data.loc[number]


'''
Count the number of intergene regions
''' 

n = 0 ; m = 0
for i in number:
    if data_new.loc[i]['INTERGENIC'] == 1:
        n += 1
    else:
        m += 1
        
        


'''
Save all SNP positions to a dictionary named data_2
''' 

        
data_1 = {}


chros = set(data_new['CHR_ID'])

for g in chros:
    if type(g) == float:
        continue
    tmp = data_new[data_new['CHR_ID'] == g]
    if ";" in g:
        keys = g.split(';')[0]
        if keys in data_1.keys():
            pass
        else:
            data_1[keys] = []
        for i in tmp.index:
            pos = tmp.loc[i]['CHR_POS']
            j = pos.split(';')
            for k in j:
                data_1[keys].append(int(k))
    elif "x" in g:
        keys = g.split(' x ')[0]
        if keys in data_1.keys():
            pass
        else:
            data_1[keys] = []
        for i in tmp.index:
            pos = tmp.loc[i]['CHR_POS']
            j = pos.split(' x ')
            for k in j:
                data_1[keys].append(int(k))
    else:
        if g in data_1.keys():
            pass
        else:
            data_1[g] = []

        for i in tmp.index:
            data_1[g].append(int(tmp.loc[i]['CHR_POS']))
 
data_2 = {}        
 
for g in data_1.keys():
    print (g , len(data_1[g]))
    d = set(data_1[g])    
    print (len(d))
    d = list(d)
    d.sort()
    data_2[g] = d


    
'''
Plot a frequency distribution histogram
'''        
        
   
        
distance = []
for g in data_2:
    a = data_2[g]
    b = [a[t+1] - a[t] for t in range(len(a) -1)]
    distance.extend(b)
    

fig = plt.figure(figsize = (12 , 10))
plt.hist(distance  , bins=50000)    
plt.xlim(0, 100000)
plt.xlabel('Fragment length', fontsize = 20)
plt.ylabel('Fragment numbers', fontsize = 20)
plt.title('Fragment length frequency distribution', fontsize = 25)
   
run_Plot(fig , 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\Fragment length frequency distribution.pdf')    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
        