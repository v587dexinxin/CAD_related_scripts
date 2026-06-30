# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:19:55 2022

@author: 86182
"""

import pandas as pd
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

def takeClosest(myList, myNumber):
    if (myNumber >= myList[-1]):
        return myList[-1]
    elif myNumber <= myList[0]:
        return myList[0]
    pos = bisect_left(myList, myNumber)   #
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before


    

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
        if ('heart' in a) or ('cardiovascular' in a) or ('coronary artery' in a) or ('cardiac' in a) or ('arrhythmia' in a) or ('hypertension' in a) or ('cardiomyopathy' in a)  or ('atherosclerosis' in a) :
            # print (i)
            number.append(i)
            break



data_new = data.loc[number]


df1=data_new.groupby(["CHR_ID","CHR_POS"]).size()
col=df1[df1>1].reset_index()[["CHR_ID","CHR_POS"]]
data_new_1 = pd.merge(col,data_new,on=["CHR_ID","CHR_POS"])
data_new_1 = data_new_1.drop_duplicates(['CHR_ID','CHR_POS','DISEASE/TRAIT'],keep='first')


data_new_2 = data_new.drop_duplicates(['CHR_ID','CHR_POS'],keep='first')



'''
Hapmap
'''

# out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD.txt' , 'w')

chro = [str(x) for x in range(1, 23)] + ['X']

SNPs = {} ; selected_SNPs = [] ; n = 0
# out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')
for g in chro:
    print (g)
    tmp_snp = pd.read_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
    tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.8]
    # SNPs[g] = tmp_snp
    tmp_data = data_new_2[data_new_2['CHR_ID'] == g]
    tmp_pos = list(tmp_data.CHR_POS)
    tmp_pos = [int(x) for x in tmp_pos]
    a = set(list(tmp_snp.pos1) + list(tmp_snp.pos2))
    a = list(a)
    a.sort()
    tmp_snp_1 = tmp_snp[tmp_snp.pos1.isin(tmp_pos) | tmp_snp.pos2.isin(tmp_pos)]
    b = set(list(tmp_snp_1.pos1) + list(tmp_snp_1.pos2))
    b = list(b)
    b.sort()
    for i in b:
        selected_SNPs.append((g , i))
        # out.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
        # n += 1
        

    for i in tmp_pos:
        if i not in a:
    #         pos = takeClosest(a, i)
    #         tmp1 = tmp_snp[tmp_snp['pos1'] == pos]
    #         tmp2 = tmp_snp[tmp_snp['pos2'] == pos]
            selected_SNPs.append((g , i))
            # out.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
            n += 1
            if i in b:
                print (i)
# out.close()
    
 
selected_SNPs = pd.DataFrame(selected_SNPs, columns=['Chromosome' , 'Position'])
selected_SNPs = selected_SNPs.sort_values(by=['Chromosome' , 'Position'])
snp_num = ['snp_' + str(x) for x in range(len(selected_SNPs))]
selected_SNPs.insert(0 , 'SNP' , snp_num , allow_duplicates=False)
selected_SNPs.to_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD0.8_all.bed' , sep = '\t' , index=False)

out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD0.8+-100bp_all.bed' , 'w')
for i in selected_SNPs.index:
    out.writelines('\t'.join(['chr' + selected_SNPs.loc[i].Chromosome , str(selected_SNPs.loc[i].Position - 100) , str(selected_SNPs.loc[i].Position + 100)]) + '\n')
out.close()


               
            
            
            
    
    
    



# selected_SNPs = [] ; n = 0
# for g in chro:
#     print (g)
#     tmp_data = data_new_2[data_new_2['CHR_ID'] == g]
#     tmp_data = list(tmp_data['SNPS'])
#     tmp_snps = SNPs[g]
#     # tmp = list(tmp_snps['rs1']) + list(tmp_snps['rs2'])
#     for i in tmp_data:
#         tmp1 = tmp_snps[tmp_snps['rs1'] == i]
#         tmp2 = tmp_snps[tmp_snps['rs2'] == i]
#         a = len(tmp1) + len(tmp2)
#         if a > 1:
#             n += 1
#         if tmp1.size > 1:
#             selected_SNPs.extend(list(tmp1.rs1))
#             selected_SNPs.extend(list(tmp1.rs2))
#         else:
#             pass
#         if tmp2.size > 1:
#             selected_SNPs.extend(list(tmp2.rs1))
#             selected_SNPs.extend(list(tmp2.rs2))     
#         else:
#             pass
        
        
        
        
    
    



















































