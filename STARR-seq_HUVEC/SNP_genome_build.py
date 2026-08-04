# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 16:49:27 2025

@author: lenovo
"""

import pandas as pd
import numpy as np
import os



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


def write_pos_to_seq_1(seqs , outfil):
    
    out = open(outfil , 'w')
    
    for i in seqs.index:
        g = seqs.loc[i]['chr']
        start = seqs.loc[i]['start'] - 1000
        end = seqs.loc[i]['end'] + 1000
        name = seqs.loc[i]['name']
        seq = hg38[g][start : end]      
        out.writelines('>' + g + '_' + str(start) + '_' + str(end) + '_' + str(name) + '\n')
            
        for j in range(len(seq) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(seq[n1:n2]).upper() + '\n')
    out.close()
        
    
    
    


hg38 = read_genome('../reference/hg38.fa')


allel = pd.read_csv('/scratch/2025-03-17/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/ref/CAD_risk_allel/CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']



risk_out = open('CAD_first6000_SNPs_risk_hg38.fa' , 'w')
nonrisk_out = open('CAD_first6000_SNPs_non_risk_hg38.fa' , 'w')

for g in chrom:
    risk_out.write('>' + g)
    nonrisk_out.write('>' + g)
    SNP = allel[allel['CHR_ID'] == g.lstrip('chr')]
    for i in range(len(hg38[g])):
        if i % 10000 == 0:
            print (g , i)
        snp_pos = i + 1
        snp = SNP[SNP['CHR_POS'] == snp_pos]
        if i % 50 == 0:
            risk_out.write('\n')
            nonrisk_out.write('\n')   
        else:
            pass
        if len(snp) > 0:
            risk = snp.iloc[0]['Risk Allel']
            nonrisk = snp.iloc[0]['Non_risk Allel']
            print (g , i , risk , nonrisk)
            risk_out.write(risk)
            nonrisk_out.write(nonrisk)
        else:
            risk_out.write(hg38[g][i])
            nonrisk_out.write(hg38[g][i])

            
            
            
risk_out.close()
nonrisk_out.close()

        
        










