# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:05:06 2023

@author: 86182
"""

import pandas as pd
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages






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
        if  ('carditis' in a) or ('heart' in a) or ('cardiovascular' in a) or ('coronary artery' in a) or ('cardiac' in a) or ('arrhythmia' in a) or ('hypertension' in a) or ('blood pressure' in a) or ('cardiomyopathy' in a)  or ('atherosclerosis' in a) or ('myocardial' in a) or ('aortic' in a):
            # print (i)
            number.append(i)
            break



data_new = data.loc[number]


df1=data_new.groupby(["CHR_ID","CHR_POS"]).size()
col=df1[df1>1].reset_index()[["CHR_ID","CHR_POS"]]
data_new_1 = pd.merge(col,data_new,on=["CHR_ID","CHR_POS"])
data_new_1 = data_new_1.drop_duplicates(['CHR_ID','CHR_POS','DISEASE/TRAIT'],keep='first')


data_new_2 = data_new.drop_duplicates(['CHR_ID','CHR_POS'],keep='first')

## replace trait ',' to '_ '
for i in data_new_2.index:
    trait = data_new_2.loc[i , 'DISEASE/TRAIT']
    if ',' in trait:
        tra = trait.replace(', ' , '_')
        print (i)
        data_new_2.loc[i , 'DISEASE/TRAIT'] = tra
        
disease_new = list(set(data_new_2['DISEASE/TRAIT']))
        
'''
Get Original disease name
'''

data1 = pd.read_table('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_name.txt' , header=0 , sep = ',' )

traits = list(set(data1['TRAIT']))


'''
Write to new files
'''


# out = open('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_name_1.csv' , 'w' , encoding='utf-8')
# out.writelines(','.join(['TRAIT' , 'Translation' , 'DISEASE']) + '\n')

df = pd.DataFrame()

for i in disease_new:
    if i in traits:
        tmp = data1[data1['TRAIT'] == i]
        
        # out.writelines(tmp.str.cat(sep=',') + '\n')
    else:
        tmp = pd.DataFrame([[i , np.nan , np.nan]])
        tmp.columns = ['TRAIT' , 'Translation' , 'DISEASE']
        
        # out.writelines(i + '\n')
    df = pd.concat([df , tmp ] , axis = 0)
        
df.to_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_name_1.csv' , index=False , sep = ','  , encoding='gbk')


    







        

