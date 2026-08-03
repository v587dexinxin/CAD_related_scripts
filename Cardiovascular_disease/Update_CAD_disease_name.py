# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:27:23 2023

@author: 86182
"""

import pandas as pd
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


data1 = pd.read_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Coronary artery disease trait_include_Hypertension.csv' , header = 0  ,  encoding='gbk')
data2 = pd.read_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Coronary artery disease trait_Hypertension_2.csv' , header = 0 )

disease_new = list(set(data2['DISEASE/TRAIT']))
traits = list(set(data1['DISEASE/TRAIT']))


'''
Write to new files
'''


# out = open('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_name_1.csv' , 'w' , encoding='utf-8')
# out.writelines(','.join(['TRAIT' , 'Translation' , 'DISEASE']) + '\n')

df = pd.DataFrame()

for i in disease_new:
    if i in traits:
        tmp = data1[data1['DISEASE/TRAIT'] == i]
        
        # out.writelines(tmp.str.cat(sep=',') + '\n')
    else:
        tmp = pd.DataFrame([[i , np.nan , np.nan]])
        tmp.columns = ['DISEASE/TRAIT', 'Translation', 'Related_TRAIT']
        
        # out.writelines(i + '\n')
    df = pd.concat([df , tmp ] , axis = 0)
        
df.to_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/Coronary artery disease trait_include_Hypertension_1.csv' , index=False , sep = ','  , encoding='gbk')


    





