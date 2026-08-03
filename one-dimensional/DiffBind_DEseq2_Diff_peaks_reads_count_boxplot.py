# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:44:32 2026

@author: lenovo
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests






chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']




deseq2 = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header = 0)
deseq2 = deseq2[deseq2['chr'].isin(chrom)]
count = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DiffBind_raw_counts_for_DESeq2.csv' , header = 0)
count = count[count['CHR'].isin(chrom)]



deseq2['LSS_count'] = 0
deseq2['OSS_count'] = 0





for g in chrom:
    tmp_deseq2 = deseq2[deseq2['chr'] == g]
    tmp_count = count[count['CHR'] == g]
    for i in tmp_deseq2.index:
        start = tmp_deseq2.loc[i]['start']
        end = tmp_deseq2.loc[i]['end']
        mask = (tmp_count['START'] <= end) & (tmp_count['END'] >= start)
        overlap = tmp_count[mask]
        if len(overlap) > 0:
             lss_count = (overlap.iloc[0]['LS_R1'] + overlap.iloc[0]['LS_R2'] + overlap.iloc[0]['LS_R3']) / 3
             oss_count = (overlap.iloc[0]['OS_R1'] + overlap.iloc[0]['OS_R2'] + overlap.iloc[0]['OS_R3']) / 3
             deseq2.loc[i , 'LSS_count'] = lss_count
             deseq2.loc[i , 'OSS_count'] = oss_count
        else:
            pass
             
                






# =========================
# 1. 读取数据
# =========================

# 如果是 csv 文件
df = deseq2.copy()

# 如果是 tab 分隔文件，改用这一行
# df = pd.read_csv("your_file.txt", sep="\t")


# =========================
# 2. 设置分组和颜色
# =========================

direction_order = [
    "LSS_specific",
    "OSS_specific",
    "not_significant"
]

direction_labels = [
    "LSS-specific",
    "OSS-specific",
    "Not significant"
]

condition_order = ["LSS", "OSS"]

colors = {
    "LSS": "#7BC96F",
    "OSS": "#C77CFF"
}


# =========================
# 3. 只保留需要展示的 peaks
# =========================

df_plot = df[df["direction"].isin(direction_order)].copy()

# log2 转换
df_plot["log2_LSS_count"] = np.log2(df_plot["LSS_count"] + 1)
df_plot["log2_OSS_count"] = np.log2(df_plot["OSS_count"] + 1)


# =========================
# 4. 转成长格式，用于画图
# =========================

df_long = df_plot.melt(
    id_vars=["peak_id", "direction"],
    value_vars=["log2_LSS_count", "log2_OSS_count"],
    var_name="condition",
    value_name="log2_normalized_count"
)

df_long["condition"] = (
    df_long["condition"]
    .str.replace("log2_", "", regex=False)
    .str.replace("_count", "", regex=False)
)


# =========================
# 5. 计算 LSS vs OSS 显著性差异
# =========================

results = []

for direction in direction_order:
    sub = df_plot[df_plot["direction"] == direction].dropna(
        subset=["log2_LSS_count", "log2_OSS_count"]
    )

    lss = sub["log2_LSS_count"]
    oss = sub["log2_OSS_count"]

    # 配对 Wilcoxon 检验
    stat, pvalue = wilcoxon(lss, oss)

    results.append({
        "direction": direction,
        "n_peaks": len(sub),
        "LSS_median": lss.median(),
        "OSS_median": oss.median(),
        "median_diff_LSS_minus_OSS": lss.median() - oss.median(),
        "pvalue": pvalue
    })

result_df = pd.DataFrame(results)

# BH 多重检验校正
result_df["padj_BH"] = multipletests(
    result_df["pvalue"],
    method="fdr_bh"
)[1]


# 显著性星号
def p_to_star(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


result_df["significance"] = result_df["padj_BH"].apply(p_to_star)

print(result_df)


# =========================
# 6. 画箱线图，自定义 LSS/OSS 间距
# =========================

plt.figure(figsize=(7, 5))
ax = plt.gca()

# 不同 peak 类型之间的距离
base_positions = np.arange(len(direction_order)) * 1.4

# LSS 和 OSS 箱线图之间的距离
offset = 0.3
# 箱体宽度
box_width = 0.5

# 记录每组最大 y 值，用于后面加显著性标注
y_max_list = []

for i, direction in enumerate(direction_order):
    group_y_values = []

    for condition in condition_order:
        data = df_long[
            (df_long["direction"] == direction) &
            (df_long["condition"] == condition)
        ]["log2_normalized_count"].dropna()

        group_y_values.extend(data.tolist())

        if condition == "LSS":
            pos = base_positions[i] - offset
        else:
            pos = base_positions[i] + offset

        bp = ax.boxplot(
            data,
            positions=[pos],
            widths=box_width,
            patch_artist=True,
            showfliers=False
        )

        for box in bp["boxes"]:
            box.set_facecolor(colors[condition])
            box.set_edgecolor("black")
            box.set_linewidth(1)

        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(1.2)

        for whisker in bp["whiskers"]:
            whisker.set_color("black")
            whisker.set_linewidth(1)

        for cap in bp["caps"]:
            cap.set_color("black")
            cap.set_linewidth(1)

    y_max_list.append(np.nanmax(group_y_values))


# =========================
# 7. 添加显著性标注
# =========================

y_range = max(y_max_list) - min(df_long["log2_normalized_count"])
line_height = y_range * 0.04
text_height = y_range * 0.06

for i, direction in enumerate(direction_order):
    x1 = base_positions[i] - offset
    x2 = base_positions[i] + offset

    y = y_max_list[i] + line_height
    star = result_df.loc[
        result_df["direction"] == direction,
        "significance"
    ].values[0]

    # 横线
    ax.plot(
        [x1, x1, x2, x2],
        [y, y + line_height, y + line_height, y],
        color="black",
        linewidth=1
    )

    # 星号
    ax.text(
        (x1 + x2) / 2,
        y + text_height,
        star,
        ha="center",
        va="bottom",
        fontsize=11
    )


# =========================
# 8. 美化图形
# =========================

ax.set_xticks(base_positions)
ax.set_xticklabels(
    direction_labels,
    rotation=25,
    ha="right"
)

ax.set_xlabel("")
ax.set_ylabel("log2(normalized count + 1)")

legend_elements = [
    Patch(facecolor=colors["LSS"], edgecolor="black", label="LSS"),
    Patch(facecolor=colors["OSS"], edgecolor="black", label="OSS")
]

ax.legend(
    handles=legend_elements,
    title="",
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\S3\\LSS_OSS_specific_peaks_boxplot_with_significance.pdf")


plt.show()



































