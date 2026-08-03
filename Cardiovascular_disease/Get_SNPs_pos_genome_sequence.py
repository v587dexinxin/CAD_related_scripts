# -*- coding: utf-8 -*-
"""
Created on Mon May 15 19:19:49 2023

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

def list_to_str(input_list):
    out_list = [x.upper() for x in input_list]
    out_list = ''.join(out_list)
    return out_list





genome_h38 = read_genome('D:/work/literature_data/genome/hg38/hg38.fa')


data = pd.read_csv('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all_risk_allel.csv' , header = 0)
data = data.sort_values(['CHR_ID' , 'CHR_POS'])

out = open('D:/work/Postdoctoral/GWAS疾病位点检测/results/CAD/first_6000/CAD_related_SNPs_LD0.99_all_risk_allel_hg38_sequence.txt' , 'w')

data_new = []

rd1 = 'CACGACGCTCTTCCGATCT'
rd2 = 'AGATCGGAAGAGCACACGT'
index1 = 'ATCACG'
index2 = 'CGATGT'

n = 0
for i in data.index:
    g = data.loc[i].CHR_ID
    pos = data.loc[i].CHR_POS
    rs = data.loc[i].SNPS
    non_risk = data.loc[i]['Non_risk Allel']
    risk = data.loc[i]['Risk Allel']
    ref = genome_h38['chr' + g][pos - 1].upper()
    
    non_seq = genome_h38['chr' + g][pos - 1 - 93 : pos - 1] + [non_risk] + genome_h38['chr' + g][pos : pos - 1 + 93]
    non_seq = rd1 + index1 + list_to_str(non_seq) + rd2
    risk_seq = genome_h38['chr' + g][pos - 1 - 93 : pos - 1] + [risk] + genome_h38['chr' + g][pos : pos - 1 + 93]
    risk_seq = rd1 + index2 + list_to_str(risk_seq) + rd2
    ref_seq = genome_h38['chr' + g][pos - 1 - 93 : pos - 1 + 93]
    ref_seq = rd1 + index1 + list_to_str(ref_seq) + rd2
    if ref != non_risk:
        n += 1
        if ref_seq[25:] != risk_seq[25:]:
            print (g , i)
    if non_risk == risk:
        print (g , i)
        break
    out.writelines(non_seq + '\n')
    out.writelines(risk_seq + '\n')
            
            
        









