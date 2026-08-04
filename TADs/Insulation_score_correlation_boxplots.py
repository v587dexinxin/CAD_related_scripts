# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 16:33:27 2026

@author: lenovo
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr




ST = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\Insulation_score\\ST_Insulation_score_25K_15.txt' , header = 0)
LSS = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\Insulation_score\\LSS_Insulation_score_25K_15.txt' , header = 0)
OSS = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\Insulation_score\\OSS_Insulation_score_25K_15.txt' , header = 0)




# 保证排序一致
ST = ST.sort_values(["chr", "Boundary"]).reset_index(drop=True)
LSS = LSS.sort_values(["chr", "Boundary"]).reset_index(drop=True)
OSS = OSS.sort_values(["chr", "Boundary"]).reset_index(drop=True)

chroms = ST["chr"].unique()

wt_lss = []
wt_oss = []
lss_oss = []

for chrom in chroms:

    st_score = ST.loc[ST["chr"] == chrom, "Insulation_Score"].values
    lss_score = LSS.loc[LSS["chr"] == chrom, "Insulation_Score"].values
    oss_score = OSS.loc[OSS["chr"] == chrom, "Insulation_Score"].values

    r1, _ = pearsonr(st_score, lss_score)
    r2, _ = pearsonr(st_score, oss_score)
    r3, _ = pearsonr(lss_score, oss_score)

    wt_lss.append(r1)
    wt_oss.append(r2)
    lss_oss.append(r3)

# 保存结果
corr_df = pd.DataFrame({
    "chr": chroms,
    "WT_vs_LSS": wt_lss,
    "WT_vs_OSS": wt_oss,
    "LSS_vs_OSS": lss_oss
})

print(corr_df)

corr_df.to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\Insulation_score\\Chromosome_InsulationScore_Correlation_15.csv",
    index=False
)





############Boxplots

# -----------------------------
# 2. 转成长表格式，便于画图
# -----------------------------
df_long = corr_df.melt(
    id_vars="chr",
    value_vars=["WT_vs_LSS", "WT_vs_OSS", "LSS_vs_OSS"],
    var_name="Comparison",
    value_name="Correlation"
)

# -----------------------------
# 3. 画箱线图 + 散点
# -----------------------------
comparisons = ["WT_vs_LSS", "WT_vs_OSS", "LSS_vs_OSS"]

plot_data = [corr_df[c].values for c in comparisons]

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
    y = corr_df[comp].values
    x = np.random.normal(loc=i, scale=0.05, size=len(y))
    ax.scatter(x, y, s=28, alpha=0.8)

ax.set_ylabel("Pearson correlation")
ax.set_xlabel("")
ax.set_title("Distribution of Insulation score correlations across chromosomes")
ax.set_ylim(0.4, 1.02)

# 去掉横向虚线/网格线
ax.grid(False)

plt.tight_layout()

plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\Insulation_score\\Chromosome_InsulationScore_Correlation_boxplot_15_1.pdf")
plt.show()

# -----------------------------
# 4. 输出每组的统计摘要
# -----------------------------
summary = corr_df[comparisons].describe().T
summary.to_csv("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\Insulation_score\\Chromosome_InsulationScore_correlation_summary.tsv", sep="\t")
print(summary)










