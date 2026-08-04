# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 11:00:27 2025

@author: lenovo
"""

import pandas as pd


def Get_enhancer_silencer(data , p , f_e , f_s):

    p_data = data[data['p.value'] <= p]
    enhancer = p_data[(p_data['Fold'] < f_e)]
    silencer = p_data[(p_data['Fold'] > f_s)]
    enhancer.columns = ['chr', 'start', 'end', 'width', 'strand', 'Conc', 'Conc_plasmid', 'Conc_cDNA', 'Fold', 'p.value', 'FDR']
    silencer.columns = ['chr', 'start', 'end', 'width', 'strand', 'Conc', 'Conc_plasmid', 'Conc_cDNA', 'Fold', 'p.value', 'FDR']
    
    return (enhancer , silencer)


def Get_SNP_names(data , SNPs):
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    data_new = []
    for g in chrom:
        tmp_SNPs = SNPs[SNPs['chr'] == g]
        tmp_data = data[data['chr'] == g]
        for i in tmp_data.index:
            s = tmp_data.loc[i]['start']
            e = tmp_data.loc[i]['end']
            mask = (tmp_SNPs['start'] <= e) & (tmp_SNPs['end'] >= s)
            overlap = tmp_SNPs[mask]
            seqname = ','.join(list(overlap['seq']))
            data_new.append(tuple(list(tmp_data.loc[i])) + tuple([seqname]))
    data_new = pd.DataFrame(data_new)
    data_new.columns = list(data.columns) + ['seq']
    
    return data_new
    
                
                
    





############Diffbind files

# data1 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HCT116_nonrisk\\HCT116_nonrisk_cDNA_VS_plasmid_deseq2.csv' , header = 0)
data1 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HCT116_nonrisk\\HCT116_nonrisk_cDNA_VS_plasmid_edgeR.csv' , header = 0)
# data2 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HCT116_risk\\HCT116_risk_cDNA_VS_plasmid_deseq2.csv' , header = 0)
data2 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HCT116_risk\\HCT116_risk_cDNA_VS_plasmid_edgeR.csv' , header = 0)
# data3 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_nonrisk\\HUVEC_nonrisk_cDNA_VS_plasmid_deseq2.csv' , header = 0)
data3 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_nonrisk\\HUVEC_nonrisk_cDNA_VS_plasmid_edgeR.csv' , header = 0)
# data4 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_risk\\HUVEC_risk_cDNA_VS_plasmid_deseq2.csv' , header = 0)
data4 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_risk\\HUVEC_risk_cDNA_VS_plasmid_edgeR.csv' , header = 0)



############SNP_files

SNPs = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\CAD_related_SNPs_LD0.99_all_risk_allel_sort.narrowPeak' , header = None , usecols = (0 , 1 , 2 , 3))
SNPs.columns = ['chr' , 'start' , 'end' , 'seq']



e1 , s1 = Get_enhancer_silencer(data1 , 0.05 , 0 , 0)
e2 , s2 = Get_enhancer_silencer(data2 , 0.05 , 0 , 0)
e3 , s3 = Get_enhancer_silencer(data3 , 0.05 , 0 , 0)
e4 , s4 = Get_enhancer_silencer(data4 , 0.05 , 0 , 0)


e1_name = Get_SNP_names(e1 , SNPs)
s1_name = Get_SNP_names(s1 , SNPs)
e2_name = Get_SNP_names(e2 , SNPs)
s2_name = Get_SNP_names(s2 , SNPs)
e3_name = Get_SNP_names(e3 , SNPs)
s3_name = Get_SNP_names(s3 , SNPs)
e4_name = Get_SNP_names(e4 , SNPs)
s4_name = Get_SNP_names(s4 , SNPs)


# e1_name.loc[e1_name["seq"].str.contains(",", na=False), "seq"]
# e2_name.loc[e2_name["seq"].str.contains(",", na=False), "seq"]
# e3_name.loc[e3_name["seq"].str.contains(",", na=False), "seq"]
# e4_name.loc[e4_name["seq"].str.contains(",", na=False), "seq"]



HCT116 = {'HCT116_nonrisk_enhancer' : e1_name , 'HCT116_nonrisk_silencer' : s1_name , 'HCT116_risk_enhancer' : e2_name , 'HCT116_risk_silencer' : s2_name}
HUVEC = {'HUVEC_nonrisk_enhancer' : e3_name , 'HUVEC_nonrisk_silencer' : s3_name , 'HUVEC_risk_enhancer' : e4_name , 'HUVEC_risk_silencer' : s4_name}



seq_HCT116_nonrisk_e = set(e1_name['seq']) - set(e2_name['seq'])
seq_HCT116_risk_e = set(e2_name['seq']) - set(e1_name['seq'])
seq_HCT116_nonrisk_s = set(s1_name['seq']) - set(s2_name['seq'])
seq_HCT116_risk_s = set(s2_name['seq']) - set(s1_name['seq'])


seq_HUVEC_nonrisk_e = set(e3_name['seq']) - set(e4_name['seq'])
seq_HUVEC_risk_e = set(e4_name['seq']) - set(e3_name['seq'])
seq_HUVEC_nonrisk_s = set(s3_name['seq']) - set(s4_name['seq'])
seq_HUVEC_risk_s = set(s4_name['seq']) - set(s3_name['seq'])




HCT116_nonrisk_e = e1_name[e1_name['seq'].isin(seq_HCT116_nonrisk_e)]
HCT116_risk_e = e2_name[e2_name['seq'].isin(seq_HCT116_risk_e)]
HCT116_nonrisk_s = s1_name[s1_name['seq'].isin(seq_HCT116_nonrisk_s)]
HCT116_risk_s = s2_name[s2_name['seq'].isin(seq_HCT116_risk_s)]


HUVEC_nonrisk_e = e3_name[e3_name['seq'].isin(seq_HUVEC_nonrisk_e)]
HUVEC_risk_e = e4_name[e4_name['seq'].isin(seq_HUVEC_risk_e)]
HUVEC_nonrisk_s = s3_name[s3_name['seq'].isin(seq_HUVEC_nonrisk_s)]
HUVEC_risk_s = s4_name[s4_name['seq'].isin(seq_HUVEC_risk_s)]









HCT116_speci = {'HCT116_nonrisk_enhancer_speci' : HCT116_nonrisk_e , 'HCT116_nonrisk_silencer_speci' : HCT116_nonrisk_s , 'HCT116_risk_enhancer_speci' : HCT116_risk_e , 'HCT116_risk_silencer_speci' : HCT116_risk_s}
HUVEC_speci = {'HUVEC_nonrisk_enhancer_speci' : HUVEC_nonrisk_e , 'HUVEC_nonrisk_silencer_speci' : HUVEC_nonrisk_s , 'HUVEC_risk_enhancer_speci' : HUVEC_risk_e , 'HUVEC_risk_silencer_speci' : HUVEC_risk_s}

speci = {'HCT116_nonrisk_enhancer_speci' : HCT116_nonrisk_e , 'HCT116_nonrisk_silencer_speci' : HCT116_nonrisk_s , 'HCT116_risk_enhancer_speci' : HCT116_risk_e , 'HCT116_risk_silencer_speci' : HCT116_risk_s , 
    'HUVEC_nonrisk_enhancer_speci' : HUVEC_nonrisk_e , 'HUVEC_nonrisk_silencer_speci' : HUVEC_nonrisk_s , 'HUVEC_risk_enhancer_speci' : HUVEC_risk_e , 'HUVEC_risk_silencer_speci' : HUVEC_risk_s}




for i in [HCT116, HUVEC]:
    for j in i:
        print (j , len(i[j]))
        i[j].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\' + j + '_p0.05_fc+-0_edgR.csv' , header = True , index = None)
        
        
    
    
for i in [HCT116_speci, HUVEC_speci]:
    for j in i:
        print (j , len(i[j]))
        i[j].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\nonrisk_VS_risk_specific_RE\\' + j + '.csv' , header = True , index = None)
        
    



for i in HCT116_speci:
    for j in HUVEC_speci:
        print (i + '_&_' + j , set(HCT116_speci[i]['seq']) & set(HUVEC_speci[j]['seq']) , len(set(HCT116_speci[i]['seq']) & set(HUVEC_speci[j]['seq'])))
        

seqs = []
seqs_dirt = {}

keys = speci.keys()
for i in keys:
    for j in keys:
        if j != i:        
            if len(set(speci[i]['seq']) & set(speci[j]['seq'])) > 0:
                seqs.extend(set(speci[i]['seq']) & set(speci[j]['seq']))
                seqs_dirt[i + '_&_' + j] = set(speci[i]['seq']) & set(speci[j]['seq'])
                print (i + '_&_' + j , set(speci[i]['seq']) & set(speci[j]['seq']) , len(set(speci[i]['seq']) & set(speci[j]['seq'])))


seqs = set(seqs)


new = []
for i in seqs:
    for j in seqs_dirt:
        if i in seqs_dirt[j]:
            new.append((i , j))
            

new = pd.DataFrame(new)



new.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\HUVEC_or_HCT116_mutative_SNPs_edgR.txt' , header = None , index = None , sep = '\t')







