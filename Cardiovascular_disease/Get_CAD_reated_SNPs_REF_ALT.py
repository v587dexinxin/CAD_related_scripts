# -*- coding: utf-8 -*-
"""
Created on Thu May 11 10:55:26 2023

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
        if len(lists) > 1 and lists[0] == '>':
            chrs = (line.split('>')[1]).split()[0]
            dict_genome[chrs] = []

        else:
            dict_genome[chrs].extend(lists)
    return dict_genome


def run_Plot(fig, OutFile):
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


genome_h38 = read_genome('H:/work/literature_data/genome/hg38/hg38.fa')
chro = [str(x) for x in range(1, 23)] + ['X']


'''
Read EBI table
'''

data = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv',
                     header=0, dtype={"REPLICATION SAMPLE SIZE": str, "CHR_ID": str, "CHR_POS": str, "SNP_ID_CURRENT": str})


'''
Select disease-related data
'''
number = []

for i in range(len(data)):
    for j in data.loc[i]:
        j = str(j)
        a = j.lower()
        if ('carditis' in a) or ('heart' in a) or ('cardiovascular' in a) or ('coronary artery' in a) or ('cardiac' in a) or ('arrhythmia' in a) or ('hypertension' in a) or ('blood pressure' in a) or ('cardiomyopathy' in a) or ('atherosclerosis' in a) or ('myocardial' in a) or ('aortic' in a):
            # print (i)
            number.append(i)
            break


data_new = data.loc[number]


df1 = data_new.groupby(["CHR_ID", "CHR_POS"]).size()
col = df1[df1 > 1].reset_index()[["CHR_ID", "CHR_POS"]]
data_new_1 = pd.merge(col, data_new, on=["CHR_ID", "CHR_POS"])
data_new_1 = data_new_1.drop_duplicates(
    ['CHR_ID', 'CHR_POS', 'DISEASE/TRAIT'], keep='first')


data_new_2 = data_new.drop_duplicates(['CHR_ID', 'CHR_POS'], keep='first')

# replace trait ',' to '_ '
for i in data_new_2.index:
    trait = data_new_2.loc[i, 'DISEASE/TRAIT']
    if ',' in trait:
        tra = trait.replace(', ', '_')
        print(i)
        data_new_2.loc[i, 'DISEASE/TRAIT'] = tra


'''
remove other disease
'''

disease = pd.read_csv(
    'H:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_name_1.csv', header=0, encoding='gbk')
all_name = set(disease.DISEASE)
disease_name = ['心血管疾病', '冠状动脉病', '高血压', '冠心病', '动脉粥样硬化', '心衰',
                '心肌病', '心肌梗死', '心脏病', '心肌肥大', '心律失常', '先天性心脏病', '心脏瓣膜病', '中风']
english_name = ['Cardiovascular disease', 'Coronary artery disease', 'Hypertension', 'Coronary heart disease',
                'Atherosclerosis', 'Heart failure', 'Cardiomyopathy', 'Myocardial infarction', 'Heart disease',
                'Myocardial hypertrophy', 'Arrhythmia', 'Congenital heart disease', 'Valvular heart disease', 'Stroke']


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

data_new_2 = data_new_2.drop(labels=index)

'''
remove unknown risk allel
'''


index = []

for i in data_new_2.index:
    allel = data_new_2.loc[i]['STRONGEST SNP-RISK ALLELE'].split('-')[1]
    if allel not in ['A', 'G', 'C', 'T']:
        index.append(i)


data_new_2 = data_new_2.drop(labels=index)


# Load CAD related SNPs (4000 tag SNP + 2000 LD 0.99 SNP)

SNPs = pd.read_table('H:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all.bed', header=0)

# db SNP annotation

vcf_snps = pd.read_csv('H:/work/Postdoctoral/GWAS疾病位点检测/results/Cardiovascular_disease/Cardiovascular_disease_related_SNPs_REF_ALT_chr.csv', header=0, dtype={'CHR_ID': str})


# Match REF, ALT and risk-allel to the selected SNPS

'''
区分6000个SNP中的ebi gwas catalog SNPs和hapmap SNPs
'''

ebi_snps = []
hapmap_snps = []
for g in chro:
    print(g)
    tmp_snps = SNPs[SNPs['CHR_ID'] == g]
    tmp_ebi_gwas = data_new_2[data_new_2['CHR_ID'] == g]
    for i in tmp_snps.index:
        pos = tmp_snps.loc[i].CHR_POS
        rs = tmp_snps.loc[i].SNPS
        overlap_ebi = tmp_ebi_gwas[tmp_ebi_gwas['CHR_POS'] == str(pos)]
        if len(overlap_ebi) != 0:
            overlap1 = overlap_ebi[overlap_ebi['SNPS'] == rs]
            if len(overlap1) != 0:
                pos1 = int(overlap1.iloc[0].CHR_POS)
                rs1 = overlap1.iloc[0].SNPS
                risk_allel = overlap1.iloc[0]['STRONGEST SNP-RISK ALLELE'].split('-')[
                    1]
                ebi_snps.append((g, pos1, rs1, '', risk_allel))
        else:
            hapmap_snps.append((g, pos, rs))

hapmap_snps = pd.DataFrame(hapmap_snps, columns=['CHR_ID', 'CHR_POS', 'SNPS'])
ebi_snps = pd.DataFrame(ebi_snps, columns=['CHR_ID', 'CHR_POS', 'SNPS', 'REF', 'ALT'])

'''
匹配 hapmap SNPs
'''

hapmap_list = []
for g in chro:
    print(g)
    tmp_snps = hapmap_snps[hapmap_snps['CHR_ID'] == g]
    tmp_vcf = vcf_snps[vcf_snps['CHR_ID'] == g]
    for i in tmp_snps.index:
        pos = tmp_snps.loc[i].CHR_POS
        rs = tmp_snps.loc[i].SNPS
        overlap_vcf = tmp_vcf[tmp_vcf['CHR_POS'] == pos]
        if len(overlap_vcf) != 0:
            overlap2 = overlap_vcf[overlap_vcf['SNPS'] == rs]
            if len(overlap2) != 0:
                pos2 = overlap2.iloc[0].CHR_POS
                rs2 = overlap2.iloc[0].SNPS
                ref = overlap2.iloc[0].REF.split(',')[0].upper()
                alt = overlap2.iloc[0].ALT.split(',')[0].upper()
                if genome_h38['chr' + g][pos2 - 1].upper() != ref:
                    print(g, i)
                hapmap_list.append((g, pos2, rs2, ref, alt))
            else:
                # print(g, i)
                pass


hapmap_list = pd.DataFrame(hapmap_list, columns=['CHR_ID', 'CHR_POS', 'SNPS', 'REF', 'ALT'])


# 检验hapmap_list中的ref序列是否正确
n = 0
for i in hapmap_list.index:
    g = hapmap_list.loc[i].CHR_ID
    pos = hapmap_list.loc[i].CHR_POS
    rs = hapmap_list.loc[i].SNPS
    ref = hapmap_list.loc[i].REF
    alt = hapmap_list.loc[i].ALT
    if genome_h38['chr' + g][pos - 1].upper() != ref.upper():
        n += 1
        print(g, i)


'''

add_snps 为data_new_classify['Hypertension'][-60:]中挑选出的risk allel 与ref 不相同SNPs。
将 add_snps中的SNPs匹配ref，并添加到需分析数据中。
'''

m = []
n = []
add_snps = ['rs11914354', 'rs5762197', 'rs6871246', 'rs7199751', 'rs7584120',
            'rs12144175', 'rs59596837', 'rs12941507', 'rs11652784', 'rs79043147', 'rs614797', 'rs258887', 'rs225748']


for i in add_snps:
    tmp = data_new_2[data_new_2.SNPS == i]
    g = tmp.iloc[0].CHR_ID
    pos = int(tmp.iloc[0].CHR_POS)
    rs = tmp.iloc[0].SNPS
    risk = tmp.iloc[0]['STRONGEST SNP-RISK ALLELE'].strip().split('-')[1].upper()
    m.append((g, pos, rs, genome_h38['chr' + g][pos - 1].upper(), risk))


'''
将ebi gwas catalog 中的 SNPs 匹配ref，并添加到需分析数据中。同时去除resk allele 与ref相同，但匹配不到ALT信息的SNPs。
'''

for g in chro:
    tmp_ebi = ebi_snps[ebi_snps['CHR_ID'] == g]
    tmp_vcf = vcf_snps[vcf_snps['CHR_ID'] == g]
    for i in tmp_ebi.index:
        pos = tmp_ebi.loc[i].CHR_POS
        rs = tmp_ebi.loc[i].SNPS
        risk = tmp_ebi.loc[i].ALT
        if genome_h38['chr' + g][pos - 1].upper() != risk.upper():
            m.append((g, pos, rs, genome_h38['chr' + g][pos - 1].upper(), risk))
        else:
            overlap = tmp_vcf[tmp_vcf['CHR_POS'] == pos]
            if len(overlap) != 0:
                overlap1 = overlap[overlap['SNPS'] == rs]
                if len(overlap1) != 0:
                    non_risk = overlap1.iloc[0].ALT.split(',')[0].upper()
                    n.append((g, pos, rs, non_risk, risk))

                else:
                    print(g, pos)
            else:
                print(g, pos)


ebi_list = pd.DataFrame(m + n, columns=['CHR_ID', 'CHR_POS', 'SNPS', 'REF', 'ALT'])


all_snps = pd.concat([ebi_list , hapmap_list] , axis = 0)
all_snps.columns = ['CHR_ID', 'CHR_POS', 'SNPS', 'Non_risk Allel', 'Risk Allel']

all_snps.to_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all_risk_allel.csv' , index = False)








