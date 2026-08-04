# -*- coding: utf-8 -*-
"""
Created on Fri May 29 20:33:21 2026

@author: lenovo
"""

import pandas as pd 
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import pandas as pd
from scipy import stats






cDNA = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\HUVEC_STARR_seq_cDNA_DESeq2_normalized_counts.txt' , header=0)
cDNA.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'HUVEC_R1_nonrisk', 'HUVEC_R2_nonrisk', 'HUVEC_R1_risk', 'HUVEC_R2_risk']


plasmid = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\HUVEC_STARR_seq_plasmid_DESeq2_normalized_counts.txt' , header=0 )
plasmid.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'HUVEC_P_R1_nonrisk', 'HUVEC_P_R2_nonrisk', 'HUVEC_P_R1_risk', 'HUVEC_P_R2_risk']


common_id = set(cDNA['Geneid']) & set(plasmid['Geneid'])
cDNA = cDNA[cDNA['Geneid'].isin(common_id)]
plasmid = plasmid[plasmid['Geneid'].isin(common_id)]

cDNA = cDNA.sort_values(by='Geneid')
plasmid = plasmid.sort_values(by='Geneid')


cDNA = cDNA.reset_index(drop=True)
plasmid = plasmid.reset_index(drop=True)



if cDNA["Geneid"].equals(plasmid["Geneid"]):
    print ('OK!')
    matrix1 = cDNA[['HUVEC_R1_nonrisk', 'HUVEC_R2_nonrisk', 'HUVEC_R1_risk', 'HUVEC_R2_risk']]
    matrix1 = np.array(matrix1.div(matrix1.sum(axis=0), axis=1))
    matrix2 = plasmid[['HUVEC_P_R1_nonrisk', 'HUVEC_P_R2_nonrisk', 'HUVEC_P_R1_risk', 'HUVEC_P_R2_risk']]
    matrix2 = np.array(matrix2.div(matrix2.sum(axis=0), axis=1))
    matrix = matrix1 / matrix2
    






matrix = pd.DataFrame(matrix)

df = pd.concat([cDNA[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']] , matrix] , axis = 1)

df.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length' , 'nonrisk_R1' , 'nonrisk_R2' , 'risk_R1' , 'risk_R2' ] 


df = df.fillna(0)

df.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\HUVEC_nonrisk_risk_cDNA_vs_plasmid_DEseq2_Norm.txt' , header = True , index = None , sep = '\t')




#################显著性计算

nonrisk_cols = ["nonrisk_R1", "nonrisk_R2"]
risk_cols = ["risk_R1", "risk_R2"]


def row_ttest(row):
    nonrisk = row[nonrisk_cols].astype(float).values
    risk = row[risk_cols].astype(float).values
    
    t_stat, p_value = stats.ttest_ind(
        nonrisk,
        risk,
        equal_var=True
    )
    
    return pd.Series({
        "t_stat": t_stat,
        "p_value": p_value
    })

df[["t_stat", "p_value"]] = df.apply(row_ttest, axis=1)

df["nonrisk_mean"] = df[nonrisk_cols].mean(axis=1)
df["risk_mean"] = df[risk_cols].mean(axis=1)

df["log2FC_risk_vs_nonrisk"] = np.log2(
    (df["risk_mean"] + 1e-6) / (df["nonrisk_mean"] + 1e-6)
)

def add_sig(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

df["significance"] = df["p_value"].apply(add_sig)

df




df.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\HUVEC_nonrisk_risk_cDNA_vs_plasmid_ttest.txt' , header = True , index = None , sep = '\t')





















