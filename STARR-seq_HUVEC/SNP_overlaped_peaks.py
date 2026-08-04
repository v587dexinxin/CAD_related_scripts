# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 19:35:31 2026

@author: lenovo
"""



import pandas as pd
import numpy as np


def SNPs_overlap_peaks(SNPs , peaks):
    SNPs['CHR_ID'] = SNPs['CHR_ID'].astype('str')
    common = []
    for g in chrom:
        tmp_snp = SNPs[SNPs['CHR_ID'] == g.lstrip('chr')]
        tmp_peaks = peaks[peaks['chr'] == g]
        for i in tmp_snp.index:
            seq = tmp_snp.loc[i]['Seq_name']
            pos = tmp_snp.loc[i]['CHR_POS']
            start = pos - 200
            end = pos + 200
            mask = (tmp_peaks['start'] <= end) & (tmp_peaks['end'] >= start)
            overlap = tmp_peaks[mask]
            if len(overlap) > 0:
                common.append((seq , g , pos))
    common = pd.DataFrame(common , columns=['Seq_name' , 'CHR_ID' , 'CHR_POS'])
    return (common)





DAVs = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\DAVs_Anno.txt' , header = 0)

common_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_common_peaks_q0.05_fc0.5_clean.bed' , header = None)
common_peaks.columns = ['chr' , 'start' , 'end']
lss_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
lss_peaks.columns = ['chr' , 'start' , 'end']
oss_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
oss_peaks.columns = ['chr' , 'start' , 'end']




common_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_OS_common_peaks_edgR_q0.05_fc1.bed' , header = None)
common_peaks.columns = ['chr' , 'start' , 'end']
lss_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_specific_peaks_edgR_q0.05_fc1.bed' , header = None)
lss_peaks.columns = ['chr' , 'start' , 'end']
oss_peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks\\HUVEC_OS_specific_peaks_deseq2.bed' , header = None)
oss_peaks.columns = ['chr' , 'start' , 'end']



deseq2_pos = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\HUVEC_LS_VS_OS_edgeR.csv' , header=0)


diff = deseq2_pos[deseq2_pos['p.value'] <= 0.05]
diff_1 = deseq2_pos[deseq2_pos['p.value'] <= 0.0001]


oss_speci_1 = diff_1[diff_1['Fold'] >= 1]
lss_speci_1 = diff[diff['Fold'] < -0.5]

no_diff = deseq2_pos[(~(deseq2_pos['p.value'] <= 0.05) ) | ((diff['Fold'] >= -0.5) & (diff['Fold'] < 1))]



common_peaks = no_diff[['seqnames' , 'start' , 'end']]
common_peaks.columns = ['chr' , 'start' , 'end']



lss_peaks = lss_speci_1[['seqnames' , 'start' , 'end']]
lss_peaks.columns = ['chr' , 'start' , 'end']


oss_peaks = oss_speci_1[['seqnames' , 'start' , 'end']]
oss_peaks.columns = ['chr' , 'start' , 'end']





chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


                

common_c = SNPs_overlap_peaks(DAVs , common_peaks)
common_l = SNPs_overlap_peaks(DAVs , lss_peaks)
common_o = SNPs_overlap_peaks(DAVs , oss_peaks)



























