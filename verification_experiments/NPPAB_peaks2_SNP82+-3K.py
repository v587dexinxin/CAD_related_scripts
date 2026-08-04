# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:14:33 2024

@author: lenovo
"""


import numpy as np
import pandas as pd
from itertools import islice



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




def write_pos_to_seq_1(POS , outfil):
    out = open(outfil , 'w')
    for i in POS:
        out.writelines('>' + i + '\n')
        for j in range(len(POS[i]) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(POS[i][n1 : n2]).upper() + '\n')
    out.close()
    
    

genome_h38 = read_genome('hg38.fa')




########### SNP82: chr1:11822085 , rs6669371 , GT

POS = {'chr1:11817085_11827085 SNP82_+-5K':genome_h38['chr1'][11817085:11827085]}

write_pos_to_seq_1(POS , '/scratch/2024-09-02/bio-shenw/Cardiovascular_disease_STARR-seq/NPPA_NPPB/chr1:11817085_11827085_SNP82_+-5K.fa')




