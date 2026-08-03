# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 17:02:31 2024

@author: lenovo
"""

import numpy as np
import pandas as pd
from itertools import islice
import os



def Load_peaks(file , peaks_type):
    sz = os.path.getsize(file)
    if sz != 0:
        peaks = pd.read_table(file , header = None)
        if peaks_type == 'narrow':
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue' , 'lengtn']
        else:
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue']
    else:
        peaks = []
    
    return peaks



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





def write_pos_to_seq(genome_h38 , g , start , end , anno , outfil):
    
    out = open(outfil , 'w')
    out.writelines('>' + g + ': ' + str(start) + '_' + str(end) + '_' + anno + '\n')
    seq = genome_h38[g][start:end]
    
    
    for j in range(len(seq) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(seq[n1 : n2]).upper() + '\n')
    out.close()
    
    




genome_h38 = read_genome('H:\\work\\literature_data\\genome\\hg38\\hg38.fa')


RNA = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)


peaks_116 = Load_peaks('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak', 'narrow')



################

PCSK9 = RNA[RNA['Gene_Name'] == 'PCSK9']

g = PCSK9['Chr'].iloc[0]
start = PCSK9['Start'].iloc[0] - 2000
end = PCSK9['Start'].iloc[0] + 1000

write_pos_to_seq(genome_h38 , g , start , end , 'PCSK9_Tss_-2K_+1K' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\PCSK9\\PCSK9_Tss_-2K_+1K.fa')


##############

peaks = peaks_116[peaks_116['name'] == 'HCT116_R2_q0.05_peak_2394']

peaks_g = peaks.iloc[0]['chr']
peaks_s = peaks.iloc[0]['start'] - 1000
peaks_e = peaks.iloc[0]['end'] + 1000


write_pos_to_seq(genome_h38 , peaks_g , peaks_s , peaks_e , 'PCSK9_related_peaks_-1K_+1K_HCT116_R2_q0.05_peak_2394' , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\PCSK9\\PCSK9_related_peaks_-1K_+1K_HCT116_R2_q0.05_peak_2394.fa')



















































