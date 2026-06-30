# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:46:43 2024

@author: lenovo
"""

import pandas as pd
import os
import numpy as np
import pandas as pd
from itertools import islice


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


def Get_SNP_peaks(SNPs , peaks):
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    n = 0; tmp = pd.DataFrame()
    for g in chrom:
        tmp_snps = SNPs[SNPs['CHR_ID'] == g.lstrip('chr')]
        tmp_peaks = peaks[peaks['chr'] == g]

        for i in tmp_snps.index:
            pos = tmp_snps.loc[i]['CHR_POS']
            
           
            overlap = tmp_peaks[(tmp_peaks['start'] <= pos + 100) & (tmp_peaks['end'] >= pos - 100)]
            if len(overlap) != 0:
                n += 1
                tmp = pd.concat([tmp , overlap])


    tmp = tmp.drop_duplicates()
    return tmp



def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['U64' , 'U64' , 'U8' , 'U4' , np.int64 , np.int64]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 0 , None):
        a = i.strip().split()
        if a[2] == 'CDS':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.split('gene_name "')[1].split('\"')[0]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf



def merge_intervals(intervals):
    if not intervals:
        return []

    # 先按照区间的开始时间排序
    intervals.sort(key=lambda x: x[0])

    merged = [intervals[0]]
    for current in intervals:
        last = merged[-1]
        if current[0] <= last[1]:  # 如果当前区间的开始时间小于等于最后一个合并区间的结束时间
            # 合并区间
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            # 不重叠的区间，直接添加到结果列表中
            merged.append(current)
    
    return merged    




SNPs = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)

peaks_huvec = Load_peaks('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional\\HUVEC_control_ls_onedimensional_q0.05_union2_peaks.narrowPeak', 'narrow')
peaks_116 = Load_peaks('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak', 'narrow')


huvec_snp_peaks = Get_SNP_peaks(SNPs , peaks_huvec)
hct116_snp_peaks = Get_SNP_peaks(SNPs , peaks_116)

union_snp_peaks = pd.concat([huvec_snp_peaks , hct116_snp_peaks])
union_snp_peaks = union_snp_peaks.drop_duplicates()
union_snp_peaks = union_snp_peaks.sort_values(by = ['chr' , 'start' , 'end'])

huvec_snp_peaks.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional\\HUVEC_SNP_peaks.narrowPeak' , header = None , index = None , sep = '\t')
hct116_snp_peaks.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional\\HCT116_SNP_peaks.narrowPeak' , header = None , index = None , sep = '\t')
union_snp_peaks.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional\\HUVEC_HCT116_union_SNP_peaks.narrowPeak' , header = None , index = None , sep = '\t')

loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\loops\\HCT116_merged6_0.1FA.hg38_loops_one_anchor_binding_loops_IF.bedpe' , header = None , usecols= (0 , 1 , 2 , 3 , 4 , 5))
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2']


RNA = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)

RNA['HCT116_FPKM'] = (RNA['HCT116_WT_R1_FPKM'] + RNA['HCT116_WT_R2_FPKM']) / 2

RNA = RNA.sort_values(by = ['HCT116_FPKM'] , ascending=False)

RNA = RNA[RNA['Chr'] != 'chrM']

gtf = Load_gtf('H:\\work\\literature_data\\genome\\hg38\\UCSC\\hg38.ncbiRefSeq.gtf\\hg38.ncbiRefSeq.gtf')


promoter = []
for i in RNA.index:
    g = RNA.loc[i]['Chr']
    strand = RNA.loc[i]['Strand']
    name = RNA.loc[i]['Gene_Name']
    fpkm = RNA.loc[i]['HCT116_FPKM']
    if strand == '+':
        start = RNA.loc[i]['Start'] - 2000
        end = RNA.loc[i]['End'] 
    elif strand == '-':
        start = RNA.loc[i]['Start']
        end = RNA.loc[i]['End'] + 2000
    else:
        print (i)
    promoter.append((g , start , end , name , fpkm))
    
promoter = pd.DataFrame(promoter)
promoter.columns = ['chr' , 'start' , 'end' , 'gene_name' , 'FPKM']




peaks_name = ['control-2_S18_q0.05_peak_347']

# selected_peak_genes = {}

for peak_name in peaks_name:
    overlap = pd.DataFrame([])
    selected_genes = pd.DataFrame([])
    peak = union_snp_peaks[union_snp_peaks['name'] == peak_name]
    g = peak.iloc[0]['chr']
    start = peak.iloc[0]['start']
    end = peak.iloc[0]['end']
    tmp_loops = loops[loops['chr1'] == g]
    mask1 = (tmp_loops['start1'] <= end) & (tmp_loops['end1'] >= start)
    mask2 = (tmp_loops['start2'] <= end) & (tmp_loops['end2'] >= start)
    overlap1 = tmp_loops[mask1]
    overlap2 = tmp_loops[mask2]
    if (len(overlap1) != 0):
        overlap = pd.concat([overlap , overlap1])
    if (len(overlap2) != 0):
        overlap = pd.concat([overlap , overlap2])    
    overlap = overlap.drop_duplicates()
    overlap.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\SNP_peaks_loops\\' + peak_name + '.bedpe' , header=None , index = None , sep = '\t')
    
 
    for i in overlap.index:
        g = overlap.loc[i]['chr1']
        start1 = overlap.loc[i]['start1']
        end1 = overlap.loc[i]['end1']
        start2 = overlap.loc[i]['start2']
        end2 = overlap.loc[i]['end2']
        tmp_promoter = promoter[promoter['chr'] == g]
        mask1 = (tmp_promoter['start'] <= end1) & (tmp_promoter['end'] >= start1)
        mask2 = (tmp_promoter['start'] <= end2) & (tmp_promoter['end'] >= start2)
        overlap1 = tmp_promoter[mask1]
        overlap2 = tmp_promoter[mask2]
        if len(overlap1) != 0:
            # print (overlap1)
            selected_genes = pd.concat([selected_genes , overlap1])
        if len(overlap2) != 0:
            # print (overlap2)
            selected_genes = pd.concat([selected_genes , overlap2])

    selected_genes = selected_genes.drop_duplicates()        
    selected_genes = selected_genes.sort_values(by = 'FPKM' , ascending=False)        
    # selected_genes = selected_genes[selected_genes['FPKM'] >= 10]
    # gene_names = list(selected_genes['gene_name'])
    # selected_peak_genes[peak_name] = gene_names
    
    
    print (selected_genes)

                




