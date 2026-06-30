# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:32:48 2022

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
        if  ('carditis' in a) or ('heart' in a) or ('cardiovascular' in a) or ('coronary artery' in a) or ('cardiac' in a) or ('arrhythmia' in a) or ('hypertension' in a) or ('cardiomyopathy' in a)  or ('atherosclerosis' in a) or ('myocardial' in a) or ('aortic' in a):
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
        tra = trait.replace(', ' , '-')
        print (i)
        data_new_2.loc[i , 'DISEASE/TRAIT'] = tra
        

'''
remove other disease
'''

disease = pd.read_table('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_disease_name.txt' , header = 0 , sep = ',')
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
disease classify
'''

data_new_classify = {'coronary artery':[] , 'hypertension':[] , 'coronary heart':[] , 'atherosclerosis':[] , 'heart failure':[] , 
                     'cardiomyopathy':[] , 'myocardial infarction':[] , 'heart disease':[] , 'myocardial hypertrophy':[] , 
                     'arrhythmia':[] , 'congenital heart disease':[] , 'valvular heart':[] , 'others':[],'stroke':[]}

# for i in data_new_2.index:
#     for j in data_new_2.loc[i]:
#         j = str(j)
#         a = j.lower()
#         for k in data_new_classify.keys():
#             if k in a:
#                 data_new_classify[k].append((data_new_2.loc[i].CHR_ID , data_new_2.loc[i].CHR_POS))
#                 break
#             else:
#                 data_new_classify[k].append((data_new_2.loc[i].CHR_ID , data_new_2.loc[i].CHR_POS))
#                 break
                






index_1 = []
data_new_classify = {}
for i in disease_name:
    data_new_classify[i] = []
    tmp = disease[disease.DISEASE == i]
    # print (i , len(tmp))
    for j in tmp.index:
        trait = tmp.loc[j].TRAIT 
        # trait = trait.replace('_' , ', ')
        tmp_1 = data_new_2[data_new_2['DISEASE/TRAIT'] == trait]
        
        # if tmp_1.size == 0:
        #     print (trait)
        for k in tmp_1.index:
            index_1.append(k)
            data_new_classify[i].append((tmp_1.loc[k].CHR_ID , tmp_1.loc[k].CHR_POS))
            
        
    
for k , v in data_new_classify.items():
    print (k , len(v))
    data_new_classify[k] = pd.DataFrame(v, columns=['CHR_ID' , 'CHR_POS'])
    
    
x = range(len(data_new_classify))
y = [len(data_new_classify[i]) for i in disease_name]
fig = plt.figure(figsize = (12 , 10))
ax = fig.add_axes([0.15  , 0.2 , 0.7 , 0.7])
plt.bar(x , y)    
plt.xlim(-0.5, 14.5)
plt.ylim(0, 5500)
plt.xlabel('Disease', fontsize = 20)
plt.xticks(x , labels = english_name , rotation = 90)
# plt.xticklabels(disease_name)
plt.ylabel('SNP numbers', fontsize = 20)
plt.title('SNP numbers associated with Cardiovascular disease', fontsize = 25)
   
run_Plot(fig , 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\cardiovascular_disease_related_SNP_Numbers_barplot_tag_SNP_new.pdf')    
   
    




'''
Hapmap
'''

# out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD.txt' , 'w')

chro = [str(x) for x in range(1, 23)] + ['X']


SNPs = {} ; selected_SNPs = {}; n = 0
# out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')

for g in chro:
    print (g)
    tmp_snp = pd.read_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
    tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.8]
    SNPs[g] = tmp_snp
        

snp_set = []
for c in disease_name:
    print (c)
    data_new_3 = data_new_classify[c]
    selected_SNPs[c] = []
    for g in chro:
        print (g)
        # tmp_snp = pd.read_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
        # tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
        tmp_snp = SNPs[g]
        # SNPs[g] = tmp_snp
        tmp_data = data_new_3[data_new_3['CHR_ID'] == g]
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
            if (g , i) not in snp_set:
                selected_SNPs[c].append((g , i))
            else:
                pass
            snp_set.append((g , i))
            # out.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
            # n += 1
            
    
        for i in tmp_pos:
            if i not in a:
        #         pos = takeClosest(a, i)
        #         tmp1 = tmp_snp[tmp_snp['pos1'] == pos]
        #         tmp2 = tmp_snp[tmp_snp['pos2'] == pos]
                selected_SNPs[c].append((g , i))
                # out.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
                n += 1
                if i in b:
                    print (i)
# out.close()
    
 
for c in selected_SNPs:
    tmp = pd.DataFrame(selected_SNPs[c], columns=['Chromosome' , 'Position'])
    tmp = tmp.sort_values(by=['Chromosome' , 'Position'])
    snp_num = ['snp_' + str(x) for x in range(len(tmp))]
    tmp.insert(0 , 'SNP' , snp_num , allow_duplicates=False)
    tmp.to_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\classify\\Cardiovascular_related_' + c + 'SNPs_LD0.8_all.bed' , sep = '\t' , index=False)
    
    out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\classify\\Cardiovascular_related_' + c + 'SNPs_LD0.8+-100bp_all.bed' , 'w')
    for i in tmp.index:
        out.writelines('\t'.join(['chr' + tmp.loc[i].Chromosome , str(tmp.loc[i].Position - 100) , str(tmp.loc[i].Position + 100)]) + '\n')
    out.close()


               
disease_name = ['心血管疾病' , '冠状动脉病' , '高血压' , '冠心病' , '动脉粥样硬化' , '心衰' , '心脏病' , '心肌病' , '心肌梗死' , '心律失常' , '中风' , '心肌肥大' , '先天性心脏病' , '心脏瓣膜病']
english_name = ['Cardiovascular disease' , 'Coronary artery disease' , 'Hypertension' , 'Coronary heart disease' ,
                'Atherosclerosis' , 'Heart failure', 'Heart disease', 'Cardiomyopathy', 'Myocardial infarction' , 
                'Arrhythmia' , 'Stroke' ,'Myocardial hypertrophy' , 'Congenital heart disease' , 'Valvular heart disease']

 
            
x = range(len(selected_SNPs))
y = [len(selected_SNPs[i]) for i in disease_name]
fig = plt.figure(figsize = (12 , 10))
ax = fig.add_axes([0.15  , 0.2 , 0.7 , 0.7])
plt.bar(x , y)    
plt.xlim(-0.5, 14.5)
plt.ylim(0, 30000)
plt.xlabel('Disease', fontsize = 20)
plt.xticks(x , labels = english_name , rotation = 90)
# plt.xticklabels(disease_name)
plt.ylabel('SNP numbers', fontsize = 20)
plt.title('SNP numbers associated with Cardiovascular disease', fontsize = 25)


for a,b in zip(x,y):   #柱子上的数字显示
 plt.text(a,b,'%d'%b,ha='center',va='bottom',fontsize=7)
 
 
 
   
run_Plot(fig , 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\cardiovascular_disease_related_SNP_Numbers_barplot_all_SNP_new.pdf')  


            
            
    
    
    



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
        
        
        
        