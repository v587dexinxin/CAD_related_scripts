# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:03:06 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



# =========================
# 1. 输入文件
# =========================
input_file = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\union_all_FPKM.csv"
out_prefix = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\RNA-seq\\FPKM\\RNAseq_FPKM_correlation"

# 如果是 csv
df = pd.read_csv(input_file)
df = df[~df["Gene_Name"].str.contains("STRG", na=False)]

# 如果是 Excel，改用：
# df = pd.read_excel("RNAseq_FPKM.xlsx")

print(df.head())
print(df.columns)

# =========================
# 2. 指定 FPKM 样本列
# 不建议完全自动 grep，避免误识别
# =========================
fpkm_cols = [
    "ST_R1_FPKM",
    "ST_R2_FPKM",
    "LS_R1_FPKM",
    "LS_R2_FPKM",
    "OS_R1_FPKM",
    "OS_R2_FPKM"
]

# 检查列是否存在
missing_cols = [c for c in fpkm_cols if c not in df.columns]
if len(missing_cols) > 0:
    raise ValueError(f"这些列在文件中没有找到: {missing_cols}")

# =========================
# 3. 提取表达矩阵
# =========================
expr = df[["Gene_Name"] + fpkm_cols].copy()

# 转为数值
for col in fpkm_cols:
    expr[col] = pd.to_numeric(expr[col], errors="coerce")

expr = expr.fillna(0)

# =========================
# 4. 去除 Gene_Name 缺失行
# =========================
expr = expr.dropna(subset=["Gene_Name"])

# =========================
# 5. 按 Gene_Name 合并重复基因
# 对 FPKM，通常用 mean 更稳妥
# =========================
expr_gene = expr.groupby("Gene_Name")[fpkm_cols].mean()

print("Before filtering:", expr_gene.shape)

# =========================
# 6. 过滤低表达基因
# 推荐至少 2 个样本 FPKM > 1
# 如果仍然相关性低，可以改成 >5
# =========================
keep = (expr_gene > 1).sum(axis=1) >= 2
expr_gene_filtered = expr_gene.loc[keep]

print("After filtering FPKM > 1 in at least 2 samples:", expr_gene_filtered.shape)

# 可选：去掉所有样本都没有变化的基因
# expr_gene_filtered = expr_gene_filtered.loc[expr_gene_filtered.var(axis=1) > 0]

print("After removing zero-variance genes:", expr_gene_filtered.shape)

# =========================
# 7. log2(FPKM + 1)
# =========================
log_expr = np.log2(expr_gene_filtered + 1)

# =========================
# 8. 重命名样本
# =========================
rename_dict = {
    "ST_R1_FPKM": "WT_Rep1",
    "ST_R2_FPKM": "WT_Rep2",
    "LS_R1_FPKM": "LS_Rep1",
    "LS_R2_FPKM": "LS_Rep2",
    "OS_R1_FPKM": "OS_Rep1",
    "OS_R2_FPKM": "OS_Rep2"
}

log_expr = log_expr.rename(columns=rename_dict)

# =========================
# 9. 计算 Pearson 相关性
# =========================
corr = log_expr.corr(method="pearson")

print("Correlation matrix:")
print(corr)

corr.to_csv(f"{out_prefix}_pearson_matrix.csv")

# =========================
# 10. 层次聚类
# =========================
distance_matrix = 1 - corr
np.fill_diagonal(distance_matrix.values, 0)

linkage_matrix = linkage(
    squareform(distance_matrix),
    method="average"
)

dendro_info = dendrogram(
    linkage_matrix,
    labels=corr.columns,
    no_plot=True
)

ordered_samples = dendro_info["ivl"]
corr_ordered = corr.loc[ordered_samples, ordered_samples]

# =========================
# 11. 画 heatmap
# =========================

def truncate_colormap(cmap_name, minval=0.0, maxval=0.88, n=256):
    cmap = plt.get_cmap(cmap_name)
    new_cmap = LinearSegmentedColormap.from_list(
        f"{cmap_name}_truncated",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

# 截断 RdYlBu，避免使用最深蓝
cmap_light_blue = truncate_colormap("RdYlBu", minval=0.0, maxval=0.9)


fig = plt.figure(figsize=(8, 8))

ax_dendro = fig.add_axes([0.08, 0.22, 0.15, 0.65])
ax_heatmap = fig.add_axes([0.25, 0.22, 0.65, 0.65])
ax_colorbar = fig.add_axes([0.25, 0.10, 0.65, 0.03])

# dendrogram
dendrogram(
    linkage_matrix,
    labels=corr.columns.tolist(),
    orientation="left",
    ax=ax_dendro,
    color_threshold=0,
    above_threshold_color="black"
)



ax_dendro.set_xticks([])
ax_dendro.set_yticks([])

for spine in ax_dendro.spines.values():
    spine.set_visible(False)

# heatmap
im = ax_heatmap.imshow(
    corr_ordered,
    aspect="auto",
    vmin=0.85,
    vmax=0.95,
    cmap=cmap_light_blue,
    origin="lower"
)

ax_heatmap.set_xticks(np.arange(len(ordered_samples)))
ax_heatmap.set_yticks(np.arange(len(ordered_samples)))

ax_heatmap.set_xticklabels(
    ordered_samples,
    rotation=45,
    ha="right",
    fontsize=11
)

ax_heatmap.set_yticklabels(
    ordered_samples,
    fontsize=11
)

ax_heatmap.yaxis.tick_right()

# grid
ax_heatmap.set_xticks(np.arange(-0.5, len(ordered_samples), 1), minor=True)
ax_heatmap.set_yticks(np.arange(-0.5, len(ordered_samples), 1), minor=True)

ax_heatmap.grid(
    which="minor",
    color="black",
    linestyle="-",
    linewidth=0.8
)

ax_heatmap.tick_params(which="minor", bottom=False, left=False)

# add numbers
for i in range(corr_ordered.shape[0]):
    for j in range(corr_ordered.shape[1]):
        value = corr_ordered.iloc[i, j]
        ax_heatmap.text(
            j,
            i,
            f"{value:.2f}",
            ha="center",
            va="center",
            color="black",
            fontsize=10
        )

# colorbar
cbar = plt.colorbar(im, cax=ax_colorbar, orientation="horizontal")
cbar.set_label("Pearson correlation", fontsize=11)

fig.suptitle(
    "Pearson Correlation of FPKM",
    fontsize=16,
    y=0.95
)

plt.savefig(f"{out_prefix}_pearson_heatmap_flipped.pdf", bbox_inches="tight")
plt.savefig(f"{out_prefix}_pearson_heatmap_flipped.png", dpi=300, bbox_inches="tight")

plt.close()