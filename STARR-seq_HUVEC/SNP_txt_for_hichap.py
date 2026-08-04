# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 00:37:37 2025

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


allel = pd.read_csv('/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/ref/CAD_risk_allel/CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']



out = open('CAD_risk_nonrisk_SNPs.txt' , 'w')

for i in allel.index:
    g = allel.loc[i]['CHR_ID']
    pos = allel.loc[i]['CHR_POS']
    risk = allel.loc[i]['Risk Allel'].upper()
    nonrisk = allel.loc[i]['Non_risk Allel'].upper()
    pos_genome = pos - 1
    ref = hg38['chr' + g][pos_genome].upper()
    out.writelines('\t'.join([g , str(pos) , ref , risk , nonrisk]) + '\n')
    
out.close()
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        