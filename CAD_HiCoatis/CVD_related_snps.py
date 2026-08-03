# -*- coding: utf-8 -*-
"""
Created on Wed May 27 11:03:48 2026

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

data['CHR_ID'] = 'chr' + data['CHR_ID']
data['start'] = data['CHR_POS'] - 200 
data['end'] = data['CHR_POS'] + 200



data[['CHR_ID' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\gwas_catalog_v1.0.2_All_diseases_SNPs+-200bp.bed' , header = None , index = None , sep = '\t')



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
        

'''
remove other disease
'''

disease = pd.read_csv('H:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_name_1.csv' , header = 0 , encoding='gbk')
all_name = set(disease.DISEASE)
disease_name = ['心血管疾病' , '冠状动脉病' , '高血压' , '冠心病' , '动脉粥样硬化' , '心衰' , '心肌病' , '心肌梗死' , '心脏病' , '心肌肥大' , '心律失常' , '先天性心脏病' , '心脏瓣膜病','中风']
english_name = ['Cardiovascular disease' , 'Coronary artery disease' , 'Hypertension' , 'Coronary heart disease' ,
                'Atherosclerosis' , 'Heart failure' , 'Cardiomyopathy', 'Myocardial infarction', 'Heart disease',
                'Myocardial hypertrophy' , 'Arrhythmia' , 'Congenital heart disease' , 'Valvular heart disease','Stroke']


other_disease_name = [x for x in all_name if x not in disease_name]
# other_disease_name.remove('心血管疾病')


index = [] 

for i in other_disease_name:
    tmp = disease[disease.DISEASE == i]
    for j in tmp.index:
        trait = tmp.loc[j].TRAIT 
        # trait = trait.replace('_' , ', ')
        tmp_1 = data_new_2[data_new_2['DISEASE/TRAIT'] == trait]
        for k in tmp_1.index:
            index.append(k)
        
data_new_2 = data_new_2.drop(labels = index)

'''
remove unknown risk allel
'''


index = []

for i in data_new_2.index:
    allel = data_new_2.loc[i]['STRONGEST SNP-RISK ALLELE'].split('-')[1]
    if allel not in ['A' , 'G' , 'C' , 'T']:
        index.append(i)


data_new_2 = data_new_2.drop(labels = index)        
    
       
    



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
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.5]
    SNPs[g] = tmp_snp
        
          
            
#### New Statistics from xxli





chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


CAD_realated_snps = []



for g in chrom:
    print(g)
    tmp_select = data_new_2[data_new_2['CHR_ID'] == g.lstrip('chr')]
    tmp_snps = SNPs[g.lstrip('chr')]
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

        CAD_realated_snps.extend(tmp)
        CAD_realated_snps.append((g , pos , rs))

                
CAD_realated_snps = pd.DataFrame(CAD_realated_snps , columns=['chr' , 'pos' , 'SNPs'])

CAD_realated_snps = CAD_realated_snps.drop_duplicates(subset=['SNPs'] , keep = 'first')


CAD_realated_snps.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CVD_related_SNPs_LD0.5_all.bed' , header = None , index = None , sep = '\t')
                        

    

    
