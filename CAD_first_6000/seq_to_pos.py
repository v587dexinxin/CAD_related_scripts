# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:53:50 2024

@author: lenovo
"""

import pandas as pd 



def seq_to_pos(element , out_fil):

    out = open(out_fil , 'w')
    for i in element.index:
        loc = int(element.loc[i]['Gene_Name'].lstrip('seq')) // 2
        g = data.loc[loc]['CHR_ID']
        pos = data.loc[loc]['CHR_POS']
        out.writelines('chr' + g + '\t' + str(pos - 1) + '\t' + str(pos + 1) + '\n')
    out.close()
    
    
    
    

data = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_sort.csv' , header = 0)


enhancer_116 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\HCT116_enhancer_deseq2_norm.csv' , header = 0 )
silencer_116 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\HCT116_silencer_deseq2_norm.csv' , header = 0 )

enhancer_huvec = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\HUVEC_enhancer_deseq2_norm.csv' , header = 0 )
silencer_huvec = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\HUVEC_silencer_deseq2_norm.csv' , header = 0 )




    
    
seq_to_pos(enhancer_116 , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\position\\HCT116_enhancer_deseq2_pos_+-1bp.bed')
seq_to_pos(silencer_116 , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\position\\HCT116_silencer_deseq2_pos_+-1bp.bed')
    
seq_to_pos(enhancer_huvec , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\position\\HUVEC_enhancer_deseq2_pos_+-1bp.bed')
seq_to_pos(silencer_huvec , 'H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\DESeq2\\enhancer_silencers\\position\\HUVEC_silencer_deseq2_pos_+-1bp.bed')




    






















