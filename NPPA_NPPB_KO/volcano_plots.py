# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 22:00:22 2025

@author: lenovo
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


diff = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\RNA_seq\\DEseq2\\HCT116_WT_VS_NP13_KO_deseq2.csv' , header = 0)


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






