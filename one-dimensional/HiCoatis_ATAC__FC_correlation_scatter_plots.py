# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 21:08:37 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import seaborn as sns





chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']




deseq2_coatis = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header = 0)
deseq2_coatis = deseq2_coatis[deseq2_coatis['chr'].isin(chrom)]


deseq2_atac = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\union_peaks_related_ATAC\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all_related_ATAC.csv' , header = 0)
deseq2_atac = deseq2_atac[deseq2_atac['chr'].isin(chrom)]





common = []
for g in chrom:
    tmp_coatis = deseq2_coatis[deseq2_coatis['chr'] == g]
    tmp_atac = deseq2_atac[deseq2_atac['chr'] == g]
    for i in tmp_coatis.index:
        start = tmp_coatis.loc[i]['start']
        end = tmp_coatis.loc[i]['end']
        direction = tmp_coatis.loc[i]['direction']
        fold = tmp_coatis.loc[i]['log2FoldChange']
        mask = (tmp_atac['start'] <= end) & (tmp_atac['end'] >= start)
        overlap = tmp_atac[mask]
        if len(overlap) > 0:
            common.append((fold , overlap.iloc[0]['log2FoldChange'] , direction))
        else:
            pass
            # common.append((fold , 0 , direction))
             
                

common = pd.DataFrame(common , columns=['Coatis_Fold' , 'ATAC_Fold' , 'direction'])



# =========================
# 2. 数据清理
# =========================
plot_df = common.copy()

plot_df["Coatis_Fold"] = pd.to_numeric(plot_df["Coatis_Fold"], errors="coerce")
plot_df["ATAC_Fold"] = pd.to_numeric(plot_df["ATAC_Fold"], errors="coerce")

plot_df = plot_df.replace([np.inf, -np.inf], np.nan)
plot_df = plot_df.dropna(subset=["Coatis_Fold", "ATAC_Fold", "direction"])

# 可选：如果只想画 differential Hi-Coatis peaks
# plot_df = plot_df[plot_df["direction"].isin(["LSS_specific", "OSS_specific"])]

# 可选：如果想画所有 peaks，包括 not_significant，就保留上面的 plot_df

# =========================
# 3. 计算相关性
# =========================
pearson_r, pearson_p = pearsonr(plot_df["Coatis_Fold"], plot_df["ATAC_Fold"])
spearman_r, spearman_p = spearmanr(plot_df["Coatis_Fold"], plot_df["ATAC_Fold"])

print(f"Pearson r = {pearson_r:.3f}, p = {pearson_p:.3e}")
print(f"Spearman r = {spearman_r:.3f}, p = {spearman_p:.3e}")

# =========================
# 4. 设置颜色
# =========================
color_dict = {
    "LSS_specific": "#6BAF75",      # 绿色
    "OSS_specific": "#D96CA6",      # 淡紫红色
    "not_significant": "#BDBDBD"    # 灰色
}

order = ["not_significant", "LSS_specific", "OSS_specific"]

# =========================
# 5. 画散点图
# =========================
plt.figure(figsize=(5.2, 4.8))

for group in order:
    sub_df = plot_df[plot_df["direction"] == group]
    
    if group == "not_significant":
        plt.scatter(
            sub_df["Coatis_Fold"],
            sub_df["ATAC_Fold"],
            s=10,
            c=color_dict[group],
            alpha=0.35,
            edgecolors="none",
            label="Not significant"
        )
    elif group == "LSS_specific":
        plt.scatter(
            sub_df["Coatis_Fold"],
            sub_df["ATAC_Fold"],
            s=16,
            c=color_dict[group],
            alpha=0.75,
            edgecolors="none",
            label="LSS-specific"
        )
    elif group == "OSS_specific":
        plt.scatter(
            sub_df["Coatis_Fold"],
            sub_df["ATAC_Fold"],
            s=16,
            c=color_dict[group],
            alpha=0.75,
            edgecolors="none",
            label="OSS-specific"
        )

# =========================
# 6. 添加整体趋势线
# =========================
sns.regplot(
    data=plot_df,
    x="Coatis_Fold",
    y="ATAC_Fold",
    scatter=False,
    color="black",
    line_kws={"linewidth": 1.5, "linestyle": "--"}
)

# =========================
# 7. 辅助线
# =========================
plt.axhline(0, color="gray", linewidth=0.8, linestyle="--")
plt.axvline(0, color="gray", linewidth=0.8, linestyle="--")

# =========================
# 8. 坐标轴和标题
# =========================
plt.xlabel("Hi-Coatis log2FC (OSS / LSS)", fontsize=12)
plt.ylabel("ATAC-seq log2FC (OSS / LSS)", fontsize=12)

plt.title(
    "Correlation between Hi-Coatis and ATAC-seq changes",
    fontsize=12
)

# =========================
# 9. 添加相关性文字
# =========================
text = (
    f"Pearson r = {pearson_r:.2f}\n"
    f"Spearman r = {spearman_r:.2f}"
)

plt.text(
    0.05,
    0.95,
    text,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top"
)

# =========================
# 10. 图例和美化
# =========================
plt.legend(
    frameon=False,
    fontsize=10,
    loc="lower right"
)

plt.tick_params(axis="both", labelsize=10)

for spine in ["top", "right"]:
    plt.gca().spines[spine].set_visible(False)

plt.tight_layout()

# =========================
# 11. 保存图片
# =========================
plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\S3\\FigS3D_HiCoatis_ATAC_correlation.pdf", dpi=300, bbox_inches="tight")

plt.show()

