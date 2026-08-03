# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 15:16:50 2024

@author: lenovo
"""

import numpy as np



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


genome_h38 = read_genome('H:/work/literature_data/genome/hg38/hg38.fa')


def write_pos_to_seq(pos , outfil):
    
    out = open(outfil , 'w')
    
    for i in pos:
        selected = genome_h38[i[0]][i[1]:i[2]]
        out.writelines('>' + i[0] + ':' + str(i[1]) + '_' + str(i[2]) + '\n')
        for j in range(len(selected) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(selected[n1 : n2]).upper() + '\n')
    out.close()
    
    
    
pos = [('chr1' , 11800000 , 11920000)]
    
write_pos_to_seq(pos , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\selected_position\\chr1_11800000_11920000.fa')  


pos = [('chr11' , 9662077 , 9664076) , ('chr11' , 102108447	, 102110446)]
write_pos_to_seq(pos , 'C:\\Users\\lenovo\\Desktop\\gaotianyu_seq.fa')  




def write_pos_to_seq_1(pos , outfil):
    
    out = open(outfil , 'w')
    
    for i in pos:
        selected = genome_h38[i[0]][i[1] - 2000 :i[2] + 2000]
        out.writelines('>' + i[0] + '_' + str(i[1]) + '_' + str(i[2]) + '_' + i[3] +  '\n')
        for j in range(len(selected) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(selected[n1 : n2]).upper() + '\n')
    out.close()    
    
pos = [('chr21' , 38171340 , 38172421 , 'KD1_seq7158') , ('chr17' , 34382438 , 34383171 , 'KD2_tar_ccl2') , ('chr4' , 73661518 , 73662220 , 'KD3_tar_cxcl8_1') , ('chr4' , 73704366 , 73705428 , 'KD4_tar_cxcl8_2')]
write_pos_to_seq_1(pos , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CRISPRI\\CRISRRI_seq.fa')  

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    