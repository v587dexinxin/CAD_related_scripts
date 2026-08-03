# -*- coding: utf-8 -*-
"""
Created on Sat May 23 14:36:28 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================
# 1. 读取差异 peak 结果
# =========================

file_path = r"H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv"

# 如果是 csv 文件
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
df = pd.read_csv(file_path)
df = df[df['chr'].isin(chrom)]

# 如果是 tab 分隔 txt 文件，用下面这一行替换上一行
# df = pd.read_csv(file_path, sep="\t")

print(df.head())
print(df.columns)


# =========================
# 2. 数据整理
# =========================

# 确保数值列是 numeric
df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
df["padj"] = pd.to_numeric(df["padj"], errors="coerce")

# 去除缺失值
df = df.dropna(subset=["log2FoldChange", "padj"])

# 避免 padj = 0 导致 -log10 无限大
df["padj_plot"] = df["padj"].replace(0, np.nextafter(0, 1))

# 计算 -log10 adjusted p-value
df["neg_log10_padj"] = -np.log10(df["padj_plot"])


# =========================
# 3. 设置阈值
# =========================

log2fc_cutoff = 0.5
padj_cutoff = 0.05

# 如果你已经有 direction 列，可以直接用
# 这里为了保险，也根据 log2FC 和 padj 重新定义一遍

df["group"] = "Not significant"

df.loc[
    (df["log2FoldChange"] >= log2fc_cutoff) & (df["padj"] < padj_cutoff),
    "group"
] = "OSS-specific"

df.loc[
    (df["log2FoldChange"] <= -log2fc_cutoff) & (df["padj"] < padj_cutoff),
    "group"
] = "LSS-specific"


# =========================
# 4. 统计数量
# =========================

print(df["group"].value_counts())

n_oss = (df["group"] == "OSS-specific").sum()
n_lss = (df["group"] == "LSS-specific").sum()
n_ns = (df["group"] == "Not significant").sum()


# =========================
# 5. 绘制火山图
# =========================

colors = {
    "LSS-specific": "#54A24B",
    "Not significant": "lightgrey",
    "OSS-specific": "#D65F9E"
}

plt.figure(figsize=(6, 5))

# 先画不显著点
for group in ["Not significant", "LSS-specific", "OSS-specific"]:
    sub_df = df[df["group"] == group]
    plt.scatter(
        sub_df["log2FoldChange"],
        sub_df["neg_log10_padj"],
        s=8,
        c=colors[group],
        alpha=0.65,
        edgecolors="none",
        label=f"{group} (n={len(sub_df)})"
    )

# 阈值线
plt.axvline(
    x=log2fc_cutoff,
    color="black",
    linestyle="--",
    linewidth=0.8
)

plt.axvline(
    x=-log2fc_cutoff,
    color="black",
    linestyle="--",
    linewidth=0.8
)

plt.axhline(
    y=-np.log10(padj_cutoff),
    color="black",
    linestyle="--",
    linewidth=0.8
)

# 图形修饰
plt.ylim([0 , 16])
plt.xlabel("log2FoldChange (OSS / LSS)", fontsize=12)
plt.ylabel("-log10(adjusted p-value)", fontsize=12)
plt.title("Differential Hi-Coatis peaks: OSS vs LSS", fontsize=14)

plt.legend(
    frameon=False,
    fontsize=9,
    loc="upper right"
)

plt.tight_layout()

# 保存
plt.savefig(
    r"H:/work/Postdoctoral/GWAS疾病位点检测/论文投稿/Figures/Fig3/HiCoatis_DiffBind_DEseq2_diff_peaks_volcano_OSS_vs_LSS.pdf",
    dpi=300,
    bbox_inches="tight"
)



plt.show()