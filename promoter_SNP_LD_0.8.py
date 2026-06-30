# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 10:12:38 2022

@author: 86182
"""

import pandas as pd
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
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


    

'''
Read promoter SNPs
'''

promoter = pd.read_table('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_related_SNPs_Promoter.txt' , header = 0 , dtype={"seqnames" : str, "start" : str , "end" : str})


'''
Hapmap
'''

# out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\Cardiovascular_related_SNPs_LD.txt' , 'w')

chro = [str(x) for x in range(1, 23)] + ['X']

SNPs = {} ; selected_SNPs = [] ; n = 0
# out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')
for g in chro:
    print (g)
    tmp_snp = pd.read_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\hapmap_hg38\\ld_chr' + g + '_CEU.bed' , header = 0 , sep = '\t')
    tmp_snp.columns = ['pos1' , 'pos2' , 'population' , 'rs1' , 'rs2' , 'Dprime' , 'R_square' , 'LOD' , 'fbin']
    tmp_snp = tmp_snp[tmp_snp['R_square'] >= 0.8]
    # SNPs[g] = tmp_snp
    tmp_data = promoter[promoter['seqnames'] == 'chr' + g]
    tmp_pos = list(tmp_data.start)
    tmp_pos = [int(x) for x in tmp_pos]
    a = set(list(tmp_snp.pos1) + list(tmp_snp.pos2))
    a = list(a)
    a.sort()
    tmp_snp_1 = tmp_snp[tmp_snp.pos1.isin(tmp_pos) | tmp_snp.pos2.isin(tmp_pos)]
    b = set(list(tmp_snp_1.pos1) + list(tmp_snp_1.pos2))
    b = list(b)
    b.sort()
    for i in b:
        selected_SNPs.append((g , i))
        # out.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
        # n += 1
        

    for i in tmp_pos:
        if i not in a:
    #         pos = takeClosest(a, i)
    #         tmp1 = tmp_snp[tmp_snp['pos1'] == pos]
    #         tmp2 = tmp_snp[tmp_snp['pos2'] == pos]
            selected_SNPs.append((g , i))
            # out.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
            n += 1
            if i in b:
                print (i)
# out.close()
    
 
selected_SNPs = pd.DataFrame(selected_SNPs, columns=['Chromosome' , 'Position'])
selected_SNPs = selected_SNPs.sort_values(by=['Chromosome' , 'Position'])
snp_num = ['snp_' + str(x) for x in range(len(selected_SNPs))]
selected_SNPs.insert(0 , 'SNP' , snp_num , allow_duplicates=False)
selected_SNPs.to_csv('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\promoter\\Cardiovascular_promoter_related_SNPs_LD0.8_all.bed' , sep = '\t' , index=False)

out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\promoter\\Cardiovascular_promoter_related_SNPs_LD0.8+-50bp_all' , 'w')
for i in selected_SNPs.index:
    out.writelines('\t'.join(['chr' + selected_SNPs.loc[i].Chromosome , str(selected_SNPs.loc[i].Position - 50) , str(selected_SNPs.loc[i].Position + 50)]) + '\n')
out.close()


out = open('D:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\Cardiovascular_disease\\promoter\\Cardiovascular_promoter_related_SNPs_LD0.8+-1bp_all' , 'w')
for i in selected_SNPs.index:
    out.writelines('\t'.join(['chr' + selected_SNPs.loc[i].Chromosome , str(selected_SNPs.loc[i].Position - 1) , str(selected_SNPs.loc[i].Position + 1)]) + '\n')
out.close()




               