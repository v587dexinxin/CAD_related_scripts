# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 10:36:53 2023

@author: 86182
"""

import pandas as pd
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



def read_genome(filename):
    file_genome = open(filename)
    dict_genome = {}
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>' :
            chrs = (line.split('>')[1]).split()[0]
            dict_genome[chrs] = []
        
        else :
            dict_genome[chrs].extend(lists)
    return dict_genome


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




genome_h38 = read_genome('D:/work/literature_data/genome/hg38/hg38.fa')


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
CAD related disease
'''

selected_disease = ['冠状动脉病' , '冠心病' , '动脉粥样硬化' , '高血压']
selected_classify = []

index_1 = []
# Data_new_classify = {}
for i in selected_disease:
    # Data_new_classify[i] = []
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
            # Data_new_classify[i].append((tmp_1.loc[k].SNPS, tmp_1.loc[k].CHR_ID , tmp_1.loc[k].CHR_POS , tmp_1.loc[k]['MAPPED_GENE'],  tmp_1.loc[k]['DISEASE/TRAIT'] , tmp_1.loc[k]['STUDY']))
            selected_classify.append((tmp_1.loc[k].SNPS, tmp_1.loc[k].CHR_ID , tmp_1.loc[k].CHR_POS , tmp_1.loc[k]['REGION'] , tmp_1.loc[k]['MAPPED_GENE'],  tmp_1.loc[k]['DISEASE/TRAIT'] , tmp_1.loc[k]['STUDY']))
    # print (i , len(Data_new_classify[i]))
    

selected_classify = pd.DataFrame(selected_classify)
selected_classify.columns = ['SNPS' , 'CHR_ID' , 'CHR_POS' , 'REGION' , 'MAPPED_GENE' , 'DISEASE/TRAIT' , 'STUDY']
t = set(selected_classify['DISEASE/TRAIT'])
t = pd.DataFrame(t)

t.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Coronary artery disease trait_Hypertension.csv' , header = ['DISEASE/TRAIT'] , index = False , sep = ',')


##注：CAD 相关trait写入之后手动进行翻译和分类！！！命名为'Coronary artery disease trait_include_Hypertension.csv'##


'''
CAD related disease classify
'''

data_new_classify = {}

CAD = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Coronary artery disease trait_include_Hypertension.csv',encoding='gbk',sep=',', header=0)


xlable = list(set(CAD.Related_TRAIT))

for i in xlable:
    data_new_classify[i] = []
    traits = CAD[CAD['Related_TRAIT'] == i]
    for j in traits.index:
        t = traits.loc[j]['DISEASE/TRAIT']
        tmp = selected_classify[selected_classify['DISEASE/TRAIT'] == t]
        for k in tmp.index:
            if not ((tmp.loc[k].CHR_ID.isdigit()) and (tmp.loc[k].CHR_ID in ['X' , 'Y'])):
                print (tmp.loc[k].CHR_ID)
                
            data_new_classify[i].append((tmp.loc[k].SNPS, tmp.loc[k].CHR_ID , tmp.loc[k].CHR_POS , tmp.loc[k]['REGION'] , tmp.loc[k]['MAPPED_GENE'],  tmp.loc[k]['DISEASE/TRAIT'] , tmp.loc[k]['STUDY']))
    data_new_classify[i] = pd.DataFrame(data_new_classify[i])
    data_new_classify[i].columns = ['SNPS' , 'CHR_ID' , 'CHR_POS' , 'REGION' , 'MAPPED_GENE' , 'DISEASE/TRAIT' , 'STUDY']


'''
Hapmap
'''

# out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD.txt' , 'w')

chro = [str(x) for x in range(1, 23)] + ['X']


SNPs = {}
# out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')

for g in chro:
    print (g)
    tmp_snp = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
    tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.99]
    SNPs[g] = tmp_snp
        


# #### New Statistics from hluo



# key_list = ['Ischemic heart disease', 'Severe coronary stenosis diseases', 'Metabolite diseases', 'Coronary artery disease in diabetes', 'Vascular endothelium diseases', 'Coronary heart disease', 'Others', 'Cholesterol lipid metabolism diseases' , 'Immune response inflammation diseases' , 'Hypertension']


# Final_Dict = {}
# UniqueSet = set()
# for key in key_list:
#     df = data_new_classify[key]
#     print(key, df.shape)
#     Final_Dict[key] = []
#     for g in chro:
#         mask = df.CHR_ID.map(lambda x : ';'+g in x or ' '+g in x or x == g)
#         sub_df = df[mask]
        
#         # 处理奇怪的pos格式
#         sub_pos = sub_df.CHR_POS.unique()
#         handled_pos = []
#         bandled_rs = []
#         for i in sub_pos:
#             if 'x' in i:
#                 i = list(map(int, i.split(' x ')))
                
#                 handled_pos.extend(i)
#             elif ';' in i:
#                 i = list(map(int, i.split(';')))
#                 handled_pos.extend(i)
#             else:
#                 i = int(i.strip())
#                 handled_pos.append(i)
#         sub_snps = SNPs[g]
        
#         mask = (sub_snps.pos1.isin(handled_pos))|(sub_snps.pos2.isin(handled_pos))
#         selected_snps = sub_snps[mask]
#         print(selected_snps.shape)
#         N = selected_snps.shape[0]
#         for i in range(N):
#             v1 = (g, selected_snps.pos1.values[i] , selected_snps.rs1.values[i])
#             v2 = (g, selected_snps.pos2.values[i] , selected_snps.rs2.values[i])
#             if v1 not in UniqueSet:
#                 Final_Dict[key].append(v1)
#                 UniqueSet.add(v1)
#             if v2 not in UniqueSet:
#                 Final_Dict[key].append(v2)
#                 UniqueSet.add(v2)
#         for i in handled_pos:
#             if (g , i) not in UniqueSet:
#                 Final_Dict[key].append((g , i))
#                 UniqueSet.add((g , i))
                
            
#### New Statistics from xxli



key_list = ['Ischemic heart disease', 'Severe coronary stenosis diseases', 'Metabolite diseases', 'Coronary artery disease in diabetes', 'Vascular endothelium diseases', 'Coronary heart disease', 'Others', 'Cholesterol lipid metabolism diseases' , 'Immune response inflammation diseases' , 'Hypertension']


Final_Dict = {}
UniqueSet = set()
for key in key_list:
    df = data_new_classify[key]
    print(key, df.shape)
    Final_Dict[key] = []
    for g in chro:
        mask = df.CHR_ID.map(lambda x : ';'+g in x or ' '+g in x or x == g)
        sub_df = df[mask]
        
        # 处理奇怪的pos格式
        sub_pos = sub_df.CHR_POS.unique()
        sub_rs = sub_df.SNPS.unique()
        if len(sub_pos) != len(sub_rs):
            print ('error')
        handled_pos = []
        for i in sub_pos:
            i = int(i.strip())
            handled_pos.append(i)
        sub_snps = SNPs[g]
        
        mask = (sub_snps.pos1.isin(handled_pos))|(sub_snps.pos2.isin(handled_pos))
        selected_snps = sub_snps[mask]
        print(selected_snps.shape)
        N = selected_snps.shape[0]
        for i in range(N):
            v1 = (g, selected_snps.pos1.values[i] , selected_snps.rs1.values[i])
            v2 = (g, selected_snps.pos2.values[i] , selected_snps.rs2.values[i])
            if v1 not in UniqueSet:
                Final_Dict[key].append(v1)
                UniqueSet.add(v1)
            if v2 not in UniqueSet:
                Final_Dict[key].append(v2)
                UniqueSet.add(v2)
        for i in range(len(handled_pos)):
            if (g , handled_pos[i] , sub_rs[i]) not in UniqueSet:
                Final_Dict[key].append((g , handled_pos[i] , sub_rs[i]))
                UniqueSet.add((g , handled_pos[i] , sub_rs[i]))
                
                        


    

    
####write to files
##classify
for c in key_list:
    tmp = pd.DataFrame(Final_Dict[c], columns=['CHR_ID' , 'CHR_POS' , 'SNPS'])
    tmp = tmp.drop_duplicates(['CHR_ID','CHR_POS'],keep='first')
    tmp = tmp.sort_values(by=['CHR_ID','CHR_POS'])
    snp_num = ['snp_' + str(x) for x in range(len(tmp))]
    tmp.insert(0 , 'SNP' , snp_num , allow_duplicates=False)
    tmp.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_' + c + 'SNPs_LD0.99_all.bed' , sep = '\t' , index=False)
    
    out = open('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_' + c + 'SNPs_LD0.99+-1bp_all.bed' , 'w')
    for i in tmp.index:
        out.writelines('\t'.join(['chr' + tmp.loc[i].CHR_ID , str(tmp.loc[i].CHR_POS - 1) , str(tmp.loc[i].CHR_POS + 1)]) + '\n')
    out.close()
       
##all    

tmp = pd.DataFrame(UniqueSet, columns=['CHR_ID' , 'CHR_POS' , 'SNPS'])
tmp = tmp.drop_duplicates(['CHR_ID','CHR_POS'],keep='first')
tmp = tmp.sort_values(by=['CHR_ID' , 'CHR_POS'])
snp_num = ['snp_' + str(x) for x in range(len(tmp))]
tmp.insert(0 , 'SNP' , snp_num , allow_duplicates=False)
tmp.to_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_SNPs_LD0.99_all.bed' , sep = '\t' , index=False)

out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_SNPs_LD0.99+-1bp_all.bed' , 'w')
for i in tmp.index:
    out.writelines('\t'.join(['chr' + tmp.loc[i].CHR_ID , str(tmp.loc[i].CHR_POS - 1) , str(tmp.loc[i].CHR_POS + 1)]) + '\n')
out.close()


out = open('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_SNPs_LD0.99+-500bp_all.bed' , 'w')
for i in tmp.index:
    out.writelines('\t'.join(['chr' + tmp.loc[i].CHR_ID , str(tmp.loc[i].CHR_POS - 500) , str(tmp.loc[i].CHR_POS + 500)]) + '\n')
out.close()


               
disease_name = ['高血压' , '糖尿病中的冠状动脉疾病' , '血管内皮疾病' , '胆固醇脂代谢疾病' , '冠心病' , '缺血性心脏病' , '严重冠状动脉狭窄疾病' ,
                '代谢疾病' , '免疫反应性炎症疾病' , '其它' ]
english_name = ['Hypertension' , 'Coronary artery disease in diabetes' , 'Vascular endothelium diseases' , 'Cholesterol lipid metabolism diseases' , 
                'Coronary heart disease' , 'Ischemic heart disease' ,  'Severe coronary stenosis diseases' , 'Metabolite diseases', 
                'Immune response inflammation diseases' , 'Others']

 
            
x = range(len(Final_Dict))
y = [len(Final_Dict[i]) for i in english_name]
fig = plt.figure(figsize = (12 , 10))
ax = fig.add_axes([0.15  , 0.2 , 0.7 , 0.7])
plt.bar(x , y)    
plt.xlim(-0.5, 9.5)
plt.ylim(0, 15000)
plt.xlabel('Disease', fontsize = 20)
plt.xticks(x , labels = english_name , rotation = 90)
# plt.xticklabels(disease_name)
plt.ylabel('SNP numbers', fontsize = 20)
plt.title('SNP numbers associated with CAD', fontsize = 25)


for a,b in zip(x,y):   #柱子上的数字显示
 plt.text(a,b,'%d'%b,ha='center',va='bottom',fontsize=7)
 
 
 
   
run_Plot(fig , 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\Plot\\CADe_related_SNP_LD0.9_Numbers_barplot_all_SNP_new_2.pdf')  



###Selected 6000 SNPs

hypertension_index = []
literature_cjin = data_new_classify['Hypertension'][data_new_classify['Hypertension'].STUDY == 'Genetic analysis of over 1 million people identifies 535 new loci associated with blood pressure traits.']
for i in literature_cjin.index:
    hypertension_index.append(i)

for i in data_new_classify['Hypertension'].index:
    if i not in hypertension_index:
        hypertension_index.append(i)
        if len(hypertension_index) == 1756:
            break
        
hypertension_index.sort()
    
data_new_classify['Hypertension'] = data_new_classify['Hypertension'].loc[hypertension_index]   




 










































