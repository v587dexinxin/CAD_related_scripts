# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 17:34:37 2026

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







def write_pos_to_seq(seq , title , outfil):
    
    out = open(outfil , 'w')
    out.writelines('>' + title + '\n')
    
        
    for j in range(len(seq) // 70 + 1):
        n1 = j * 70
        n2 = (j + 1) * 70
        out.writelines(''.join(seq[n1:n2]).upper() + '\n')
    out.close()


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
        
        



hg38 = read_genome('/scratch/2026-05-02/bio-shenw/ref/Human/hg38/hg38.fa')



#####9p21
seq = hg38['chr9'][21950000 : 22140000]
seq1 = ''.join(seq)

write_pos_to_seq(seq1 , 'chr9:21950000-22140000_CAD_risk_region' , '/scratch/2026-04-13/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/chr9:21950000-22140000_CAD_risk_region.fasta')  





#####
seq = hg38['chr17'][34980000 : 34990000]
seq1 = ''.join(seq)

write_pos_to_seq(seq1 , 'chr17:34980000-34990000_seq4758_CCL2_related_SNPs' , '/scratch/2026-04-13/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/chr17:34980000-34990000_seq4758_CCL2_related_SNPs.fasta')  






























