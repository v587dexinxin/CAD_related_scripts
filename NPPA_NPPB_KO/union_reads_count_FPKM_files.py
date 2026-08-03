# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:21:57 2025

@author: lenovo
"""



import numpy as np
import pandas as pd
import os


#######HCT116_NPPA_NPPB_KO


'''
Load files
'''

data_count = {}

keys = os.listdir('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\reads_count\\')


for k in keys:
    new_k = k.split('-')[1]
    file = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\reads_count\\' + k  , sep='\t' , index_col = 0 ,header = 0 , skiprows=1)
    file = file[~file.index.duplicated(keep='first')]
    file.columns = ['Chr' , 'Start' , 'End' , 'Strand' , 'Length' , new_k + '_Count' ]
    data_count[new_k] = file

    
data_FPKM = {}
keys = os.listdir('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\FPKM\\')
new_keys = []
for k in keys:
    new_k = k.split('-')[1]
    new_keys.append(new_k)
    file = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\FPKM\\' + k , sep='\t' , index_col = 0 ,header = 0 )
    file = file[~file.index.duplicated(keep='first')]
    file.columns = ['Gene_ID' , 'Chr' , 'Strand' , 'Start' , 'End' , 'Coverage' , new_k + '_FPKM' , 'TPM']
    data_FPKM[new_k] = file
    
# data_FPKM['WT-1']['Gene_ID'] = data_FPKM['WT-1'].index

union = pd.concat([data_FPKM[new_keys[0]]['Chr'] , data_FPKM[new_keys[0]]['Strand'] ,
                   data_FPKM[new_keys[0]]['Start'] , data_FPKM[new_keys[0]]['End'] , 
                   data_count[new_keys[0]][new_keys[0] + '_Count'] , data_count[new_keys[1]][new_keys[1] + '_Count'] , 
                   data_count[new_keys[2]][new_keys[2] + '_Count'] , data_count[new_keys[3]][new_keys[3] + '_Count'] , 
                   data_count[new_keys[4]][new_keys[4] + '_Count'] , data_count[new_keys[5]][new_keys[5] + '_Count'] , 
                   data_count[new_keys[6]][new_keys[6] + '_Count'] , data_count[new_keys[7]][new_keys[7] + '_Count'] , 
                   data_count[new_keys[8]][new_keys[8] + '_Count'] , data_count[new_keys[9]][new_keys[9] + '_Count'] , 
                   data_count[new_keys[10]][new_keys[10] + '_Count'] , data_count[new_keys[11]][new_keys[11] + '_Count'] , 
                   data_count[new_keys[12]][new_keys[12] + '_Count'] , data_count[new_keys[13]][new_keys[13] + '_Count']] , axis = 1 , join = 'inner')


union.to_csv(os.path.join('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\reads_count\\NPPA_NPPB_peak2_KO_union_all_reads_count.csv') , index_label = 'Gene_Name' , sep = ',')


union = pd.concat([data_FPKM[new_keys[0]]['Chr'] , data_FPKM[new_keys[0]]['Strand'] ,
                   data_FPKM[new_keys[0]]['Start'] , data_FPKM[new_keys[0]]['End'] , 
                   data_FPKM[new_keys[0]][new_keys[0] + '_FPKM'] , data_FPKM[new_keys[1]][new_keys[1] + '_FPKM'] , 
                   data_FPKM[new_keys[2]][new_keys[2] + '_FPKM'] , data_FPKM[new_keys[3]][new_keys[3] + '_FPKM'] , 
                   data_FPKM[new_keys[4]][new_keys[4] + '_FPKM'] , data_FPKM[new_keys[5]][new_keys[5] + '_FPKM'] ,
                   data_FPKM[new_keys[6]][new_keys[6] + '_FPKM'] , data_FPKM[new_keys[7]][new_keys[7] + '_FPKM'] ,
                   data_FPKM[new_keys[8]][new_keys[8] + '_FPKM'] , data_FPKM[new_keys[9]][new_keys[9] + '_FPKM'] ,
                   data_FPKM[new_keys[10]][new_keys[10] + '_FPKM'] , data_FPKM[new_keys[11]][new_keys[11] + '_FPKM'] ,
                   data_FPKM[new_keys[12]][new_keys[12] + '_FPKM'] , data_FPKM[new_keys[13]][new_keys[13] + '_FPKM']] , axis = 1 , join = 'inner')


union.to_csv(os.path.join('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\FPKM\\NPPA_NPPB_peak2_KO_union_all_FPKM.csv') , index_label = 'Gene_Name' , sep = ',')







