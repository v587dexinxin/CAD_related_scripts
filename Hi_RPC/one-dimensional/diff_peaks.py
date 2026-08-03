# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 15:16:42 2025

@author: lenovo
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

diff_deseq2 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\HUVEC_LS_VS_OS_deseq2.csv' , header = 0)
diff_edgeR = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\HUVEC_LS_VS_OS_edgeR.csv' , header = 0)

diff_deseq2 = diff_deseq2[diff_deseq2['seqnames'].isin(chrom)]
diff_edgeR = diff_edgeR[diff_edgeR['seqnames'].isin(chrom)]


diff = diff_edgeR[diff_edgeR['FDR'] <= 0.05]
diff_ls = diff[diff['Fold'] <= -1]
diff_os = diff[diff['Fold'] >= 1]

common1 = diff[(diff['Fold'] > -1) & (diff['Fold'] < 1)]
common2 = diff_edgeR[diff_edgeR['FDR'] > 0.05]
common = pd.concat([common1 , common2])



diff_ls.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_specific_peaks_edgR_q0.05_fc1.csv' , header = True , index = None)
diff_os.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_OS_specific_peaks_edgR_q0.05_fc1.csv' , header = True , index = None)
common.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_OS_common_peaks_edgR_q0.05_fc1.csv' , header = True , index = None)



diff_ls[['seqnames' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_specific_peaks_edgR_q0.05_fc1.bed' , header = None , index = None , sep = '\t')
diff_os[['seqnames' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_OS_specific_peaks_edgR_q0.05_fc1.bed' , header = None , index = None , sep = '\t')
common[['seqnames' , 'start' , 'end']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\HUVEC_LS_OS_common_peaks_edgR_q0.05_fc1.bed' , header = None , index = None , sep = '\t')






# 设置阈值
fc_thresh = 1

p_thresh = 0.05

# 计算 -log10(p-value)
diff["neg_log10_FDR"] = -np.log10(diff["FDR"])

# 定义显著性
diff["significant"] = (
    (diff["FDR"] < p_thresh) & (abs(diff["Fold"]) > fc_thresh)
)

# 画图
plt.figure(figsize=(6, 6))

# 非显著点
plt.scatter(
    diff.loc[~diff["significant"], "Fold"],
    diff.loc[~diff["significant"], "neg_log10_FDR"],
    s=10,
    alpha=0.5,
    label="Not significant"
)

# 显著点
plt.scatter(
    diff.loc[diff["significant"], "Fold"],
    diff.loc[diff["significant"], "neg_log10_FDR"],
    s=15,
    alpha=0.8,
    label="Significant"
)

# 阈值线
plt.axhline(-np.log10(p_thresh), linestyle="--")
plt.axvline(fc_thresh, linestyle="--")
plt.axvline(-fc_thresh, linestyle="--")

# 标签
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(FDR)")
plt.title("Volcano Plot")

plt.legend(frameon=False)
plt.tight_layout()
plt.show()


















