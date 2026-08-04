# -*- coding: utf-8 -*-
"""
Created on Sun May 31 17:29:00 2026

@author: lenovo
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\HUVEC_nonrisk_risk_cDNA_vs_plasmid_ttest.txt' , header = 0 , sep = '\t')
DAVs = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\DAVs_Anno.txt' , header = 0)
DAVs = list(DAVs['Seq_name'])


# 确保数值列是 numeric
df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
df["log2FC_risk_vs_nonrisk"] = pd.to_numeric(
    df["log2FC_risk_vs_nonrisk"],
    errors="coerce"
)

# 去掉缺失值
plot_df = df.dropna(subset=["p_value", "log2FC_risk_vs_nonrisk"]).copy()

# 防止 p_value = 0 导致 -log10 无穷大
plot_df["p_value_for_plot"] = plot_df["p_value"].replace(0, 1e-300)

# 计算 -log10(p-value)
plot_df["minus_log10_p"] = -np.log10(plot_df["p_value_for_plot"])

# 判断是否为 DAVs
plot_df["is_DAV"] = plot_df["Geneid"].isin(DAVs)

# 统计 DAV 数目
dav_count = plot_df["is_DAV"].sum()

# 开始画图
plt.figure(figsize=(8, 6))

# 非 DAVs，灰色
plt.scatter(
    plot_df.loc[~plot_df["is_DAV"], "log2FC_risk_vs_nonrisk"],
    plot_df.loc[~plot_df["is_DAV"], "minus_log10_p"],
    s=28,
    c="gray",
    alpha=0.3,
    edgecolors="none"
)

# DAVs，粉红色
plt.scatter(
    plot_df.loc[plot_df["is_DAV"], "log2FC_risk_vs_nonrisk"],
    plot_df.loc[plot_df["is_DAV"], "minus_log10_p"],
    s=35,
    c="lightcoral",
    alpha=0.9,
    edgecolors="none"
)

# p = 0.05 阈值线，可选
plt.axhline(
    -np.log10(0.05),
    color="gray",
    linestyle="--",
    linewidth=1,
    alpha=0.6
)

# x = 0 参考线，可选
plt.axvline(
    0,
    color="gray",
    linestyle="--",
    linewidth=1,
    alpha=0.6
)

# 标题，类似示例图
plt.text(
    0.5,
    0.98,
    f"DAVs:{dav_count}",
    transform=plt.gca().transAxes,
    ha="center",
    va="top",
    fontsize=16
)

plt.xlabel("log2(fold change)", fontsize=16)
plt.ylabel("-log10(p-value)", fontsize=16)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# 根据你的数据范围可调整
plt.xlim(-2.5, 2.5)
plt.ylim(-0.5, 6)

# 去掉上边和右边框
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\DAVs_volcano_plot.pdf", dpi=300)
plt.show()