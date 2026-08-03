# -*- coding: utf-8 -*-
"""
Created on Thu May 11 15:43:09 2023

@author: 86182
"""

import pandas as pd
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



chro = [str(x) for x in range(1, 23)] + ['X']


'''
Read EBI table
'''

data = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv' , header = 0 , dtype={"REPLICATION SAMPLE SIZE" : str, "CHR_ID" : str , "CHR_POS" : str , "SNP_ID_CURRENT" : str})

data = data[['CHR_ID' , 'CHR_POS' , 'SNPS']]    
data = data[data['CHR_ID'].isin(chro)]
data['CHR_POS'] = data['CHR_POS'].astype(int)
data = data.drop_duplicates(keep = 'first')


'''
Hapmap data
'''

SNPs = {}
# out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')

for g in chro:
    print (g)
    tmp_snp = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
    tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.9]
    SNPs[g] = tmp_snp
        
    
hapmap_snps = pd.DataFrame(columns = ['CHR_ID' , 'CHR_POS' , 'SNPS'])
    
for k , v in SNPs.items():
    print (k)
    v1 = v.drop_duplicates(['pos1' , 'rs1'] , keep = 'first')
    v2 = v.drop_duplicates(['pos2' , 'rs2'] , keep = 'first')
    
    chr_id = [k for x in range(len(v1) + len(v2))]
    chr_id = pd.Series(chr_id)
    chr_pos = pd.concat([v1.pos1 , v2.pos2])
    chr_pos.index = range(0 , len(chr_pos))
    chr_snps = pd.concat([v1.rs1 , v2.rs2])
    chr_snps.index = range(0 , len(chr_snps))
    df = pd.concat([chr_id , chr_pos , chr_snps] , axis = 1)
    df = df.drop_duplicates(keep = 'first')
    df.columns = ['CHR_ID' , 'CHR_POS' , 'SNPS']
    
    hapmap_snps = pd.concat([hapmap_snps , df] , axis = 0)
            
'''
All SNPS to be matched to VCF
'''    
snps_tobe_matched = pd.concat([data , hapmap_snps] , axis = 0)    
            
snps_tobe_matched = snps_tobe_matched.sort_values(by = ['CHR_ID' , 'CHR_POS'])            
snps_tobe_matched = snps_tobe_matched.drop_duplicates(keep = 'first')        
        
        
        
###检验结果    
# n = 0    
# for g in chro:
#     print (g)
#     tmp1 = data[data['CHR_ID'] == g].values
#     # tmp2 = SNPs[g]
#     tmp3 = snps_tobe_matched[snps_tobe_matched['CHR_ID'] == g].values     
#     for i in tmp1:
#         n += 1
#         if n % 10000 == 0:
#             print (n)
#         if i not in tmp3:
#             print (g , i)
            
        
        
        
            
# n = 0            
# for g in chro:
#     print (g)
#     # tmp1 = data[data['CHR_ID'] == g]
#     tmp2 = SNPs[g]
#     tmp3 = snps_tobe_matched[snps_tobe_matched['CHR_ID'] == g].values        
#     for i in tmp2.index:
#         n += 1
#         if n % 10000 == 0:
#             print (n)
#         pos1 = tmp2.loc[i].pos1
#         pos2 = tmp2.loc[i].pos2
#         rs1 = tmp2.loc[i].rs1
#         rs2 = tmp2.loc[i].rs1
#         if ([g , pos1 , rs1] not in tmp3) or ([g , pos2 , rs2] not in tmp3):
#             print (g , i)
                   
        
        
'''
VCF data
'''   

# n = 0  
# selected_vcf = []
# for g in chro:
#     vcf = open('D:/work/literature_data/genome/hg38/SNP/GCF_000001405.40_ChrNC_000001.11.vcf' , 'r')
#     tmp_snps = data[data['CHR_ID'] == g]
#     for i in vcf:
#         n += 1
#         if n %10000 == 0:
#             print (n)
#         i = i.strip().split()
#         c = i[0].lstrip('NC_00000').split('.')[0]
#         try:
#             pos = int(i[1])
#         except:
#             pos = i[1]
#         snps = i[2]
#         ref = i[3]
#         alt = i[4]
#         tmp = [c , pos , snps , ref , alt]
#         overlap = tmp_snps[(tmp_snps['CHR_POS'] == pos)]
#         if len(overlap) > 0:
#             selected_vcf.append(tmp)
 
            
selected_vcf = pd.DataFrame(columns = ['CHR_ID' , 'CHR_POS' , 'SNPS' , 'REF' , 'ALT'])
for g in chro:   
    print (g)         
    vcf = pd.read_table('H:/work/literature_data/genome/hg38/SNP/GCF_000001405.40_ChrNC_0000' + g + '.vcf' , header=None)
    tmp_snps = snps_tobe_matched[snps_tobe_matched['CHR_ID'] == g]
    vcf.columns = ['CHR_ID' , 'CHR_POS' , 'SNPS' , 'REF' , 'ALT']
    vcf.index = vcf.CHR_POS
    b = vcf.index
    a = list(tmp_snps.CHR_POS)
    c = []
    for i in a:
        if i in b:
            #a.remove(i)
            c.append(i)

    #         x = i -1
    #         if x in vcf.index:
    #             a.append(x)
    # for i in a:
    #     if i not in vcf.index:
    #         a.remove(i)
    selected_vcf = pd.concat([selected_vcf , vcf.loc[c]])
    

selected_vcf.to_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_related_SNPs_REF_ALT.csv' , index = False)            
vcf_chrs = [x.lstrip('NC_00000').split('.')[0] for x in selected_vcf.CHR_ID]

for i in range(len(vcf_chrs)):
    if vcf_chrs[i] == '23':
        vcf_chrs[i] = 'X'
    else:
        pass

selected_vcf.CHR_ID = vcf_chrs
selected_vcf.to_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_related_SNPs_REF_ALT_chr.csv' , index = False)   
    
    
    
    
    

            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        