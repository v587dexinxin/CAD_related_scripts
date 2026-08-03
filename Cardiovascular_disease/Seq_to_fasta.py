# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 17:37:02 2023

@author: 86182
"""
import pandas as pd 

pos_data = pd.read_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all_risk_allel_sort.csv' , header = 0)
data = open('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all_risk_allel_hg38_sequence.txt' , 'r')
out = open('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all_risk_allel_hg38_sequence.fa' , 'w')
rd1 = 'CACGACGCTCTTCCGATCT'
rd2 = 'AGATCGGAAGAGCACACGT'
index1 = 'ATCACG'
index2 = 'CGATGT'

n = 0
for i in data:
    i = i[19:211]
    index = n // 2
    chro = pos_data.loc[index].CHR_ID
    pos = pos_data.loc[index].CHR_POS
    non_risk = pos_data.loc[index]['Non_risk Allel']
    risk = pos_data.loc[index]['Risk Allel']
    out.writelines('>seq' + str(n) + ' chr' + chro + ':' + str(pos - 93) + '-' +  str(pos + 92) + '_' + non_risk + risk + '\n')
    out.writelines(i + '\n')
    n += 1
out.close()
    
    
    