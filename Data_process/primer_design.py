# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 20:51:21 2024

@author: lenovo
"""

import pandas as pd
import numpy as np


def GC_percent(seq):
    a = 0 ; g = 0 ; c = 0 ; t = 0
    for i in seq:
        if i == 'A':
            a += 1
        elif i == 'G':
            g += 1
        elif i == 'C':
            c += 1
        else:
            t += 1
    p = np.round((g + c) / (a + g + c + t) , 2)
    
    print (p)
    return p
    
def reverse_complement(dna_sequence):
    # 将输入的DNA序列转换为大写，以便处理小写字母
    dna_sequence = dna_sequence.upper()

    # 构建互补字典
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # 生成反向互补序列
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reversed(dna_sequence))

    return reverse_complement_sequence

        
    
    

data = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\selected_silencer_enhancer.csv' , header = 0)

seq_fa = open('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_hg38_sequence.fa' , 'r')

seq_data = {}
for i in seq_fa:
    if '>' in i:
        s_name = i.split()[0].lstrip('>')
    else:
        seq_data[s_name] = i.strip()
        
seq_fa.close()


selected_data={}
for i in data.columns:
    selected_data[i] = {}
    for j in data[i]:
        # s_name = 'seq' + str(int(j.lstrip('seq')) + 1)
        selected_data[i][j] = seq_data[j]
        # selected_data[i][s_name] = seq_data[s_name]
        
        
kpnI = 'GGTACC'
mluI = 'ACGCGT'


selected_primer={}
for i in selected_data:
    prefix = i.split('_')[0][0] + i.split('_')[1][0]
    # print (prefix)
    for j in selected_data[i]:
        keys = prefix + '_' + j.lstrip('seq')
        print (keys)
        f = kpnI + selected_data[i][j][:20]
        r = mluI + reverse_complement(selected_data[i][j][172:])
        f_p = GC_percent(f)
        r_p = GC_percent(r)
        if f_p < 0.4:
            f = 'GCC' + f
        elif f_p > 0.6:
            f = 'AAT' + f
        else:
            f = 'AAC' + f
        if r_p < 0.4:
            r = 'GCC' + r
        elif r_p > 0.6:
            r = 'ATT' + r
        else:
            r = 'AAC' + r
        GC_percent(f)
        GC_percent(r)
        selected_primer[keys + '_F'] = f
        selected_primer[keys + '_R'] = r
        
selected_primer = pd.DataFrame(list(selected_primer.items()), columns=['Primer', 'Sequence'])


selected_primer.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\HCT116_HUVEC_enhancer_silencer_primer_design.csv' , header = True , index=None)
        
        
        

        
                                      
                               
        
        
        
        
        
        
        
        