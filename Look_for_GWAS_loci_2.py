# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 10:13:56 2022

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

data = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv' , header = 0 , dtype={"REPLICATION SAMPLE SIZE" : str, "CHR_ID" : str , "CHR_POS" : str , "SNP_ID_CURRENT" : str})


'''
Select disease-related data
'''
number = []

for i in range(len(data)):
    for j in data.loc[i]:
        j = str(j)
        a = j.lower()
        if ('carditis' in a) or ('heart' in a) or ('cardiovascular' in a) or ('coronary artery' in a) or ('cardiac' in a) or ('arrhythmia' in a) or ('hypertension' in a) or ('cardiomyopathy' in a)  or ('atherosclerosis' in a) or ('myocardial' in a) or ('aortic' in a):
            # print (i)
            number.append(i)
            break



data_new = data.loc[number]


##Get all of SNPs duplicated
df1=data_new.groupby(["CHR_ID","CHR_POS"]).size()
col=df1[df1>1].reset_index()[["CHR_ID","CHR_POS"]]
data_new_1 = pd.merge(col,data_new,on=["CHR_ID","CHR_POS"])
##Obtain SNPs with unique characteristics
data_new_1 = data_new_1.drop_duplicates(['CHR_ID','CHR_POS','DISEASE/TRAIT'],keep='first')

##Get all uniq SNPs
data_new_2 = data_new.drop_duplicates(['CHR_ID','CHR_POS'],keep='first')




'''
Count the number of intergene regions
''' 

n = 0 ; m = 0
for i in data_new_2.index:
    if data_new_2.loc[i]['INTERGENIC'] == 1:
        n += 1
    elif data_new_2.loc[i]['INTERGENIC'] == 0:
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

keys = [str(t) for t in range(1,23)] + ['X']
# keys.sort()


outfil = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease_related_SNP_+-1Kb.bed' , 'w')
for g in keys:
    for i in data_2[g]:
        outfil.writelines('chr' + g + '\t' + str(i - 1000) + '\t' + str(i + 1000) + '\n')
outfil.close()

outfil = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease_related_SNP_+-1bp.bed' , 'w')
for g in keys:
    for i in data_2[g]:
        outfil.writelines('chr' + g + '\t' + str(i - 1) + '\t' + str(i + 1) + '\n')
outfil.close()



outfil = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease_related_SNP.bed' , 'w')
outfil.writelines('\t'.join(['SNP','Chromosome','Position']) + '\n')
n = 0
for g in keys:
    for i in data_2[g]:
        outfil.writelines('snp_' + str(n) + '\t'  + g + '\t' + str(i) + '\n')
        n += 1
outfil.close()



distance = []
data_3 = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease_related_SNP_+-1Kb_merged.bed' , 'r')
for i in data_3.readlines():
    i = i.split('\t')
    distance.append(int(i[2]) - int(i[1]))
data_3.close()
    
        
'''
Plot a frequency distribution histogram
'''        
        
   
        
# distance = []
# for g in data_2:
#     a = data_2[g]
#     b = [a[t+1] - a[t] for t in range(len(a) -1)]
#     distance.extend(b)
    

fig = plt.figure(figsize = (12 , 10))
plt.hist(distance  , bins=1000)    
plt.xlim(1900, 4000)
# plt.ylim(0, 300)
plt.xlabel('Fragment length', fontsize = 20)
plt.ylabel('Fragment numbers', fontsize = 20)
plt.title('Fragment length frequency distribution', fontsize = 25)
   
run_Plot(fig , 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\Tag_SNP_Fragment length frequency distribution.pdf')    
   
    
   



