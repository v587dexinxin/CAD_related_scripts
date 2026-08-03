# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 17:58:50 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
from itertools import islice




#########loops#############
loops_wt = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops_wt.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']
loops_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops_ls.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']
loops_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bedpe' , header = None)
loops_os.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']


# loops_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops_chr4.bedpe' , header = None)
# loops_os.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']




###########selected_activate_regulatory_elements

selected_31_51 = open('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\Selected_31-51.fasta' , 'r')

for i in selected_31_51:
    if '>' in i:
        ids = i.split('_')
        # s = ids[0].lstrip('>')
        g = ids[1]
        start = int(ids[2]) - 1000
        end = int(ids[3]) + 1000
        n = i.lstrip('>').rstrip('\n')

        tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
        tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
        tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


        overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
        overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
        overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





        overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\WT_' + n + '.bedpe' , header = None , index = None , sep = '\t')
        overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\LS_' + n + '.bedpe' , header = None , index = None , sep = '\t')
        overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\OS_' + n + '.bedpe' , header = None , index = None , sep = '\t')


        
        


selected_31_51.close()





#########loops_IF > 2#############

loops_wt = loops_wt[loops_wt['IF'] > 2]
loops_ls = loops_ls[loops_ls['IF'] > 2]
loops_os = loops_os[loops_os['IF'] > 2]

###########selected_activate_regulatory_elements

selected_31_51 = open('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\Selected_31-51.fasta' , 'r')

for i in selected_31_51:
    if '>' in i:
        ids = i.split('_')
        # s = ids[0].lstrip('>')
        g = ids[1]
        start = int(ids[2]) - 1000
        end = int(ids[3]) + 1000
        n = i.lstrip('>').rstrip('\n')

        tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
        tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
        tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


        overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
        overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
        overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





        overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\IF_2+\\WT_' + n + '.bedpe' , header = None , index = None , sep = '\t')
        overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\IF_2+\\LS_' + n + '.bedpe' , header = None , index = None , sep = '\t')
        overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_activate_enhancer_silencer_31_51\\selected_loops\\IF_2+\\OS_' + n + '.bedpe' , header = None , index = None , sep = '\t')


        
        


selected_31_51.close()






###########selected_CCL2

genes = [('chr17' , 34255285 , 34257203 , 'CCL2') , ('chr8' , 125430358 , 125438403 , 'TRIB1') , ('chr9' , 21967752 , 21995324 , 'CDKN2A') , ('chr9' , 22002903 , 22009313 , 'CDKN2B')]

for i in genes:
    g = i[0]
    start = i[1] - 2000
    end = i[2] + 2000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_WT.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_LS.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_OS.bedpe' , header = None , index = None , sep = '\t')


    
    


###########selected_1

SNPs = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)

selected = ['seq7158' , 'seq9454']
selected = ['seq82']

for i in selected:
    
    snp = SNPs[SNPs['Seq_name'] == i]
    g = 'chr' + snp.iloc[0]['CHR_ID']
    start = snp.iloc[0]['CHR_POS'] - 1000
    end = snp.iloc[0]['CHR_POS'] + 1000
    

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_loops\\WT_' + i + '.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_loops\\LS_' + i + '.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_loops\\OS_' + i + '.bedpe' , header = None , index = None , sep = '\t')


    
    







###########selected_CXCL8

genes = [('chr4' , 73740569 , 73743716 , 'CXCL8') , ('chr7' , 101127104 , 101139247 , 'SERPINE1') , ('chr6' , 12256484 , 12297194 , 'EDN1') , ('chr17' , 34382465 - 2000 , 34382864 + 2000 , 'ccl2_s')]
# genes = [('chr4' , 73740569 , 73743716 , 'CXCL8')]


for i in genes:
    g = i[0]
    start = i[1] - 2000
    end = i[2] + 2000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_WT.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_LS.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_OS_1.bedpe' , header = None , index = None , sep = '\t')


    
###########selected_EDN1

genes = [('chr6' , 12256484 , 12297194 , 'EDN1')]

for i in genes:
    g = i[0]
    start = i[1] - 2000
    end = i[2] + 2000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_WT.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_LS.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_OS.bedpe' , header = None , index = None , sep = '\t')


    

###########selected_KLF4, KLF2

genes = [('chr9' , 107484852 , 107489769 , 'KLF4')] 
# genes = [('chr19' , 16324826 , 16328685 , 'KLF2')]
genes = [('chr1' , 11907725 , 11909104 , 'NP_SE')]


for i in genes:
    g = i[0]
    start = i[1] - 5000
    end = i[2] + 5000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_WT.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_LS.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_OS.bedpe' , header = None , index = None , sep = '\t')


    



###########selected_NP_Peak1

genes = [('chr1' , 11805659 , 11806450 , 'NP_Peak1') , ('chr1' , 11821986 , 11822185 , 'seq82')] 


union1 = pd.DataFrame([]) ; union2 = pd.DataFrame([]) ; union3 = pd.DataFrame([])
for i in genes:
    g = i[0]
    start = i[1] - 1000
    end = i[2] + 1000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]
    
    union1 = pd.concat([union1 , overlap1])
    union2 = pd.concat([union2 , overlap2])
    union3 = pd.concat([union3 , overlap3])





union1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_seq82_WT.bedpe' , header = None , index = None , sep = '\t')
union2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_seq82_LS.bedpe' , header = None , index = None , sep = '\t')
union3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_seq82_OS.bedpe' , header = None , index = None , sep = '\t')




###########selected_NP_Peak1_peak2_peak3

genes = [('chr1' , 11805659 , 11806450 , 'NP_Peak1') , ('chr1' , 11821986 , 11822185 , 'seq82') , ('chr1' , 11907725 , 11909104 , 'NP_Peak3')] 


union1 = pd.DataFrame([]) ; union2 = pd.DataFrame([]) ; union3 = pd.DataFrame([])
for i in genes:
    g = i[0]
    start = i[1] - 1000
    end = i[2] + 1000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]
    
    union1 = pd.concat([union1 , overlap1])
    union2 = pd.concat([union2 , overlap2])
    union3 = pd.concat([union3 , overlap3])





union1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_peak2_peak3_WT.bedpe' , header = None , index = None , sep = '\t')
union2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_peak2_peak3_LS.bedpe' , header = None , index = None , sep = '\t')
union3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_peak2_pesk3_OS.bedpe' , header = None , index = None , sep = '\t')


    




    

###########selected_NP_Peak1_HCT116

genes = [('chr1' , 11805659 , 11806450 , 'NP_Peak1') , ('chr1' , 11821986 , 11822185 , 'seq82')] 
# genes = [('chr19' , 16324826 , 16328685 , 'KLF2')]
loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\loops\\HCT116_merged6_0.1FA.hg38_loops_one_anchor_binding_union_peaks.bedpe' , header = None)
loops.columns = ['chr1' , 's1' , 'e1' , 'chr2' , 's2' , 'e2' , 'IF' , 'qvalue']


union1 = pd.DataFrame([]) ; union2 = pd.DataFrame([]) ; union3 = pd.DataFrame([])
for i in genes:
    g = i[0]
    start = i[1] - 1000
    end = i[2] + 1000
    gene_name = i[3]

    tmp_wt = loops[(loops['chr1'] == g) & (loops['chr2'] == g)]



    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]

    
    union1 = pd.concat([union1 , overlap1])


union1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_NP_Peak1_seq82_HCT116_WT.bedpe' , header = None , index = None , sep = '\t')





###########selected_DUSP5

genes = [('chr10' , 110497907 , 110511533 , 'DUSP5')]

for i in genes:
    g = i[0]
    start = i[1] - 2000
    end = i[2] + 2000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_WT.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_LS.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_OS.bedpe' , header = None , index = None , sep = '\t')



###########selected_SELE

genes = [('chr1' , 169722640 , 169734079 , 'SELE')]

for i in genes:
    g = i[0]
    start = i[1] - 2000
    end = i[2] + 2000
    gene_name = i[3]

    tmp_wt = loops_wt[(loops_wt['chr1'] == g) & (loops_wt['chr2'] == g)]
    tmp_ls = loops_ls[(loops_ls['chr1'] == g) & (loops_ls['chr2'] == g)]
    tmp_os = loops_os[(loops_os['chr1'] == g) & (loops_os['chr2'] == g)]


    overlap1 = tmp_wt[((tmp_wt['s1'] <= end) & (tmp_wt['e1'] >= start)) | ((tmp_wt['s2'] <= end) & (tmp_wt['e2'] >= start))]
    overlap2 = tmp_ls[((tmp_ls['s1'] <= end) & (tmp_ls['e1'] >= start)) | ((tmp_ls['s2'] <= end) & (tmp_ls['e2'] >= start))]
    overlap3 = tmp_os[((tmp_os['s1'] <= end) & (tmp_os['e1'] >= start)) | ((tmp_os['s2'] <= end) & (tmp_os['e2'] >= start))]





    overlap1.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_WT.bedpe' , header = None , index = None , sep = '\t')
    overlap2.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_LS.bedpe' , header = None , index = None , sep = '\t')
    overlap3.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\selected_genes\\selected_' + gene_name + '_OS.bedpe' , header = None , index = None , sep = '\t')
