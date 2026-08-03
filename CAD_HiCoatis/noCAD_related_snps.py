# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:14:54 2026

@author: lenovo
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

data = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv' , header = 0 , dtype={"REPLICATION SAMPLE SIZE" : str, "CHR_ID" : str , "CHR_POS" : str , "SNP_ID_CURRENT" : str})
data = data[pd.to_numeric(data['CHR_POS'], errors="coerce").notna()]
data = data.dropna(subset=['CHR_POS'])
data['CHR_POS'] = data['CHR_POS'].astype(int)
data = data.drop_duplicates(subset=['SNPS'] , keep = 'first')


data['CHR_ID'] = 'chr' + data['CHR_ID']
data['start'] = data['CHR_POS'] - 200 
data['end'] = data['CHR_POS'] + 200



CAD = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_SNPs_LD0.8_all.bed' , header = None)
CAD.columns=['chr' , 'pos' , 'SNPs']
data = data[~data['SNPS'].isin(list(CAD['SNPs']))]


data[['CHR_ID' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\gwas_catalog_v1.0.2_noCAD_diseases_SNPs+-200bp.bed' , header = None , index = None , sep = '\t')


'''
Hapmap (选择LD大于0.99的SNP)
'''

# out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD.txt' , 'w')

chro = [str(x) for x in range(1, 23)] + ['X']


SNPs = {}
# out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')

for g in chro:
    print (g)
    tmp_snp = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
    tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.8]
    SNPs[g] = tmp_snp
        
          
            
#### New Statistics from xxli





chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']






for g in chrom:
    print(g)
    tmp_select = data[data['CHR_ID'] == g]
    tmp_snps = SNPs[g.lstrip('chr')]
    tmp_1 = []
    CAD_realated_snps = []
    for i in tmp_select.index:
        rs = tmp_select.loc[i]['SNPS']
        pos = tmp_select.loc[i]['CHR_POS']
        mask = (tmp_snps['rs1'] == rs) | (tmp_snps['rs2'] == rs) 
        overlap = tmp_snps[mask]
        if len(overlap) > 0:
            tmp = []
            for j in overlap.index:
                rs1 = overlap.loc[j]['rs1']
                pos1 = overlap.loc[j]['pos1']
                rs2 = overlap.loc[j]['rs2']
                pos2 = overlap.loc[j]['pos2']
                tmp.append((g , pos1 , rs1))
                tmp.append((g , pos2 , rs2))
                tmp = list(set(tmp))

        tmp_1.extend(tmp)
        tmp_1.append((g , pos , rs))
    tmp_1 = list(set(tmp_1))
    CAD_realated_snps.extend(tmp_1)
    

                
    CAD_realated_snps = pd.DataFrame(CAD_realated_snps , columns=['chr' , 'pos' , 'SNPs'])
    
    CAD_realated_snps = CAD_realated_snps.drop_duplicates(subset=['SNPs'] , keep = 'first')
    
    
    CAD_realated_snps.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\tmp\\CVD_related_SNPs_LD0.8_all_' + g + '.bed' , header = None , index = None , sep = '\t')
                            
    
        
noCAD = pd.DataFrame([] , columns = ['chr' , 'pos' , 'SNPs'])
for g in chrom:
    noCAD_chr = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\tmp\\CVD_related_SNPs_LD0.8_all_' + g + '.bed' , header = None)
    noCAD_chr.columns = ['chr' , 'pos' , 'SNPs']
    noCAD = pd.concat([noCAD , noCAD_chr])
    
    
    
noCAD = noCAD.drop_duplicates(subset=['SNPs'] , keep = 'first')    
    
    
    
    
noCAD.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\noCAD_related_SNPs_LD0.8_all_background.bed' , header = None , index = None , sep = '\t')


noCAD['start'] = noCAD['pos'] - 200
noCAD['end'] = noCAD['pos'] + 200
    
    
noCAD[['chr' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\noCAD_related_SNPs_LD0.8_all_background_+-200bp.bed' , header = None , index = None , sep = '\t')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    