# -*- coding: utf-8 -*-
"""
Created on Wed May  6 21:15:14 2026

@author: lenovo
"""

import matplotlib.pyplot as plt
import pandas as pd




loops_wt = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bedpe')
loops_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe')
loops_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bedpe')



# =========================
# 1. 输入 loops 数目
# =========================
# 请将下面的数值替换成你自己的 loops 数目
data = {
    "Condition": ["WT", "LSS", "OSS"],
    "Loops": [len(loops_wt), len(loops_ls), len(loops_os)]
}

df = pd.DataFrame(data)

# =========================
# 2. 设置颜色
# =========================
colors = {
    "WT": "#4C78A8",   # 蓝色
    "LSS": "#59A14F",  # 绿色
    "OSS": "#E15759"   # 红色
}

bar_colors = [colors[x] for x in df["Condition"]]

# =========================
# 3. 绘制柱状图
# =========================
plt.figure(figsize=(4.5, 4))

bars = plt.bar(
    df["Condition"],
    df["Loops"],
    color=bar_colors,
    width=0.6,
    edgecolor="black",
    linewidth=0.8
)

# =========================
# 4. 添加数值标签
# =========================
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        height,
        f"{int(height):,}",
        ha="center",
        va="bottom",
        fontsize=11
    )

# =========================
# 5. 美化图形
# =========================
plt.ylabel("Number of loops", fontsize=12)
plt.xlabel("")
plt.title("Hi-Coatis loops under different shear stress conditions", fontsize=12)

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.spines["top"].set_visible(False)
plt.spines["right"].set_visible(False)

plt.tight_layout()

# =========================
# 6. 保存图片
# =========================
plt.savefig("WT_LSS_OSS_loop_counts_barplot.pdf", dpi=300)
plt.savefig("WT_LSS_OSS_loop_counts_barplot.png", dpi=300)

plt.show()