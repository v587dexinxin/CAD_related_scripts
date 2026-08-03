# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:06:45 2026

@author: lenovo
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 读取 E1 矩阵
df = pd.read_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\WT_LS_OS_500k_E1_matrix.tsv", sep="\t")

# 去除 nan
df = df.replace(["nan", "NaN", "NA", "."], pd.NA)
df[["WT", "LS", "OS"]] = df[["WT", "LS", "OS"]].apply(pd.to_numeric, errors="coerce")
df = df.dropna(subset=["WT", "LS", "OS"])

#手动检查方向性
df.loc[df["chrom"] == "chr15", "OS"] *= -1
df.loc[df["chrom"] == "chr17", "LS"] *= -1
df.loc[df["chrom"] == "chr20", "OS"] *= -1
df.loc[df["chrom"] == "chr21", "WT"] *= -1
df.loc[df["chrom"] == "chr22", "LS"] *= -1
df.to_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\WT_LS_OS_500k_E1_matrix_double_checking.tsv.tsv", sep="\t")


# Pearson 相关性：适合比较 E1 数值整体一致性
pearson_corr = df[["WT", "LS", "OS"]].corr(method="pearson")
print("Pearson correlation:")
print(pearson_corr)

# Spearman 相关性：适合比较排序/趋势一致性，受极端值影响较小
spearman_corr = df[["WT", "LS", "OS"]].corr(method="spearman")
print("Spearman correlation:")
print(spearman_corr)

# 保存结果
pearson_corr.to_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\WT_LS_OS_E1_Pearson_correlation.tsv", sep="\t")
spearman_corr.to_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\WT_LS_OS_E1_Spearman_correlation.tsv", sep="\t")

# 画 Pearson heatmap
plt.figure(figsize=(4, 3.5))
sns.heatmap(
    pearson_corr,
    annot=True,
    fmt=".3f",
    cmap="RdBu_r",
    vmin=-1,
    vmax=1,
    square=True,
    linewidths=0.5,
    cbar_kws={"label": "Pearson r"}
)
plt.title("Correlation of compartment E1")
plt.tight_layout()
plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\WT_LS_OS_E1_Pearson_correlation_heatmap.pdf")
# plt.savefig("WT_LS_OS_E1_Pearson_correlation_heatmap.png", dpi=300)
plt.close()

# 画 Spearman heatmap
plt.figure(figsize=(4, 3.5))
sns.heatmap(
    spearman_corr,
    annot=True,
    fmt=".3f",
    cmap="RdBu_r",
    vmin=-1,
    vmax=1,
    square=True,
    linewidths=0.5,
    cbar_kws={"label": "Spearman rho"}
)
plt.title("Correlation of compartment E1")
plt.tight_layout()
plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\WT_LS_OS_E1_Spearman_correlation_heatmap.pdf")
# plt.savefig("WT_LS_OS_E1_Spearman_correlation_heatmap.png", dpi=300)
plt.close()











###############boxplots
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

data = []

for g in chrom:
    tmp = df[df['chrom'] == g]
    pearson_corr = tmp[["WT", "LS", "OS"]].corr(method="pearson")
    print(g, "Pearson correlation:")
    print(pearson_corr)
    data.append([g , pearson_corr['WT']['LS'] , pearson_corr['WT']['OS'] , pearson_corr['LS']['OS']])


df = pd.DataFrame(data, columns=["chrom", "WT_vs_LS", "WT_vs_OS", "LS_vs_OS"])




# 保存原始表
df.to_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\chromosome_compartment_correlations.tsv", sep="\t", index=False)

# -----------------------------
# 2. 转成长表格式，便于画图
# -----------------------------
df_long = df.melt(
    id_vars="chrom",
    value_vars=["WT_vs_LS", "WT_vs_OS", "LS_vs_OS"],
    var_name="Comparison",
    value_name="Correlation"
)

# -----------------------------
# 3. 画箱线图 + 散点
# -----------------------------
comparisons = ["WT_vs_LS", "WT_vs_OS", "LS_vs_OS"]

plot_data = [df[c].values for c in comparisons]

fig, ax = plt.subplots(figsize=(6, 5))

bp = ax.boxplot(
    plot_data,
    labels=comparisons,
    patch_artist=True,
    widths=0.55,
    showfliers=False
)

for patch in bp["boxes"]:
    patch.set_alpha(0.6)

for median in bp["medians"]:
    median.set_linewidth(1.5)

np.random.seed(123)
for i, comp in enumerate(comparisons, start=1):
    y = df[comp].values
    x = np.random.normal(loc=i, scale=0.05, size=len(y))
    ax.scatter(x, y, s=28, alpha=0.8)

ax.set_ylabel("Pearson correlation")
ax.set_xlabel("")
ax.set_title("Distribution of compartment correlations across chromosomes")
ax.set_ylim(0.6, 1.02)

# 去掉横向虚线/网格线
ax.grid(False)

plt.tight_layout()
# plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\chromosome_compartment_correlation_boxplot.png", dpi=300)
plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\chromosome_compartment_correlation_boxplot.pdf")
plt.show()

# -----------------------------
# 4. 输出每组的统计摘要
# -----------------------------
summary = df[comparisons].describe().T
summary.to_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\chromosome_compartment_correlation_summary.tsv", sep="\t")
print(summary)





