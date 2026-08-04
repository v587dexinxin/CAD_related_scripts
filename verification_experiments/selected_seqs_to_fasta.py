# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 14:52:57 2025

@author: lenovo
"""

import pandas as pd



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




def write_pos_to_seq(peaks , flanking , outfil):
    
    out = open(outfil , 'w')
    
    for i in peaks:
        g = i[0]
        start = i[1]
        end = i[2]
        f_start = start - flanking
        f_end = end + flanking
        seq = genome_h38[g][f_start:f_end]
        classify = i[4]
        peak_name = i[6]
        out.writelines('>' + g + ': ' + str(start) + '_' + str(end) + '_' +  peak_name + '_flanking+-' + str(flanking) + '_' + classify + '\n')
        for j in range(len(seq) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(seq[n1 : n2]).upper() + '\n')
    out.close()
    
    
def wirte_selecetd_seqs_to_fasta(selected_seqs , flanking , out_fil):
    
    out = open(out_fil , 'w')
    n = 1
    for i in selected_seqs:
        s = SNPs[SNPs['Seq_name'] == i]
        g = 'chr' + s.iloc[0]['CHR_ID']
        pos = s.iloc[0]['CHR_POS']
        start = pos - flanking
        end = pos + flanking
        snps = s.iloc[0]['SNPS']
        g_seq = genome_h38[g][start:end]
        if i == 'seq5074':
            i_1 = 'seq5074,seq5076,seq5078'
        else:
            i_1 = i
        p = mutative_SNPs[mutative_SNPs['seqs'] == i_1].iloc[0]['property']
        out.writelines('>s' + str(n) + '_' + g + '_' + str(pos) + '_+-' + str(flanking // 1000) + 'kb_' + i_1 + '_' + snps + '_' + p + '\n')
        out.writelines(''.join(g_seq).upper() + '\n\n')
        n += 1
        
    out.close()    


def wirte_selecetd_seq_peaks_to_fasta(selected_seqs , act_regulatory , flanking , cells , out_fil):
    
    out = open(out_fil , 'w')
    n = 31
    for i in selected_seqs:
        
        if ',' in i:
            i_1 = i.split(',')[0]
        else:
            i_1 = i
        s = SNPs[SNPs['Seq_name'] == i_1]
        act_e = act_regulatory[act_regulatory['seq'] == i]
        g = act_e.iloc[0]['chr']
        pos = s.iloc[0]['CHR_POS']
        start = act_e.iloc[0]['peaks_start'] - flanking
        end = act_e.iloc[0]['peaks_end'] + flanking
        snps = s.iloc[0]['SNPS']
        g_seq = genome_h38[g][start:end]

        if act_e.iloc[0]['nonrisk_e'] == 1:
            p1 = 'nonrisk_enhancer'
        else:
            p1 = ''
        if act_e.iloc[0]['nonrisk_s'] == 1:
            p2 = 'nonrisk_silencer'
        else:
            p2 = ''
        if act_e.iloc[0]['risk_e'] == 1:
            p3 = 'risk_enhancer'
        else:
            p3 = ''
        if act_e.iloc[0]['risk_s'] == 1:
            p4 = 'risk_silencer'
        else:
            p4 = ''  
        p = p1 + p2 + p3 + p4
        out.writelines('>s' + str(n) + '_' + g + '_' + str(start) + '_' + str(end) + '_peaks_+-' + str(flanking // 1000) + 'kb_' \
                       + i + '_' + snps + '_' + str(pos) + '_' + cells + '_' + p + '_speci\n')
        out.writelines(''.join(g_seq).upper() + '\n\n')
        n += 1
        
    out.close()    



    



genome_h38 = read_genome('/scratch/2026-01-12/bio-shenw/ref/Human/hg38/hg38.fa')

RNA = pd.read_csv('/scratch/2025-10-27/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/RNA_seq/mapping/bam2/DEGs/union_all_FPKM.csv' , header = 0)
SNPs = pd.read_csv('/scratch/2025-10-27/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)

mutative_SNPs = pd.read_table('/scratch/2025-10-27/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/HUVEC_or_HCT116_mutative_SNPs.txt' , header = None)
mutative_SNPs.columns = ['seqs' , 'property']


####selected_seqs
'''risk 和 nonrisk发生变化，并且在HCT116和HUVEC两个细胞中均有变化的seqs，一共53个'''

selected_seqs = ['seq6260', 'seq852', 'seq960', 'seq10260', 'seq1204', 'seq1620', 'seq3546', 'seq4800', 'seq4864', 'seq6442', 'seq6472', 'seq7650', 
          'seq8490', 'seq9108', 'seq942', 'seq10266', 'seq1070', 'seq10940', 'seq290', 'seq404', 'seq4340', 'seq5272', 'seq5920', 'seq6676', 
          'seq7432', 'seq836', 'seq8510', 'seq10382', 'seq2544', 'seq6998', 'seq9152', 'seq3594', 'seq318', 'seq6214', 'seq536', 'seq10422',
          'seq1962', 'seq5300', 'seq2540', 'seq1216', 'seq4808', 'seq4592', 'seq4726', 'seq2340', 'seq10776', 'seq8586', 'seq11204', 'seq3954',
          'seq10452', 'seq6482', 'seq5074', 'seq1386', 'seq9996']

###selected_1
'''HUVEC_nonrisk 或者 HCT116_nonrisk 特异性变化优先，带HCT116 Hi-Coatis一维信号的seqs优先验证'''
selected_1 = ['seq6260', 'seq852', 'seq960', 'seq10260', 'seq1204', 'seq1620', 'seq3546', 'seq4800', 'seq4864', 'seq6442', 'seq6472', 'seq7650',
              'seq8490', 'seq9108', 'seq942', 'seq5920', 'seq7432', 'seq8510', 'seq10382', 'seq2544', 'seq6998', 'seq9152', 'seq3594', 'seq318',
              'seq4592', 'seq2340', 'seq3954', 'seq1386', 'seq9996' , 'seq5074']


selected_2 = []

for i in selected_seqs:
    if i not in selected_1:
        selected_2.append(i)
        print (i)



out = '/scratch/2025-10-27/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/Selected_regulatory_element_verification_all.fasta'
out1 = '/scratch/2025-10-27/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/Selected_regulatory_element_verification_1.fasta'
out2 = '/scratch/2025-10-27/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/Selected_regulatory_element_verification_2.fasta'





wirte_selecetd_seqs_to_fasta(selected_seqs , 1000 , out)
wirte_selecetd_seqs_to_fasta(selected_1 , 1000 , out1)
wirte_selecetd_seqs_to_fasta(selected_2 , 1000 , out2)





##########selected_HUVEC_activate_enhancer

HUVEC_act_enhancer = pd.read_csv('/scratch/2026-02-02/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/activate_regulatory_elements/HUVEC_activate_nonrisk_specific_enhancer.csv' , header = 0)
HUVEC_act_silencer = pd.read_csv('/scratch/2026-02-02/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/activate_regulatory_elements/HUVEC_activate_nonrisk_specific_silencer.csv' , header = 0)
SNPs = pd.read_csv('/scratch/2026-02-02/bio-shenw/Cardiovascular_disease_STARR-seq/SNP_related_fragments/CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)


HUVEC_act_regulatory = pd.concat([HUVEC_act_enhancer , HUVEC_act_silencer])

selectd_HUVEC_act_regulatory = list(HUVEC_act_enhancer['seq']) + list(HUVEC_act_silencer['seq'])

selectd_HUVEC_act_regulatory.remove('seq4800')
selectd_HUVEC_act_regulatory.remove('seq2544')


out = '/scratch/2026-01-12/bio-shenw/Cardiovascular_disease_STARR-seq/verification_experiments/Selected_HUVEC_activate_regulatory_element_verification31-51.fasta'

wirte_selecetd_seq_peaks_to_fasta(selectd_HUVEC_act_regulatory , HUVEC_act_regulatory , 1000 , 'HUVEC' , out)

