# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 16:11:26 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



FRiP = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\FRiP\\HUVEC_Hi-Coatis_WT_LSS_OSS_FRiP.txt' , header=0)
FRiP = FRiP[FRiP['Types'] == 'Inside']

peaks_wt = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\peaks\\union_peaks\\HiRPC_WT_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2 , 4) , header = None)
peaks_wt.columns = ['chr' , 'start' , 'end' , 'score']

peaks_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\peaks\\union_peaks\\HiRPC_LS_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks_ls.columns = ['chr' , 'start' , 'end']


peaks_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\peaks\\union_peaks\\HiRPC_OS_allreps_q0.05_peaks_sorted_merged.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks_os.columns = ['chr' , 'start' , 'end']





# 数据
methods = list(FRiP['labels'])
frip_values = list(FRiP['Percentage'])
peaks = [len(peaks_coatis), len(peaks_pro) , len(peaks_rchip) , len(peaks_polII)]



# 颜色和大小
colors = plt.cm.Blues(np.interp(frip_values, (min(frip_values), max(frip_values)), (0.2, 1)))
sizes = [frip * 20 for frip in frip_values]  # 大小与frip_values成正比

# 创建图形
fig, ax = plt.subplots(figsize=(8, 5))

# 绘制散点图
for i in [3 , 2 , 1 , 0]:
    method = methods[i]
    ax.hlines(y=method, xmin=0, xmax=peaks[i], color='black', linestyles='solid', linewidth=1)
    ax.scatter(peaks[i], methods[i], s=sizes[i], color=colors[i], label=f"FRiP: {frip_values[i]:.2f}", edgecolors="black", zorder=2)
    
    
# # 添加颜色条
# sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin=min(frip_values), vmax=max(frip_values)))
# sm.set_array([])
# cb = fig.colorbar(sm, ax=ax, orientation='horizontal', label='FRiP Values', shrink=0.6, pad=0.2)



# 添加图例和标题
legend = ax.legend(title="Methods", loc="lower right", bbox_to_anchor=(1.55, 0), handletextpad=0.5, borderpad=0.5, labelspacing=1.5)
ax.set_title("Peak Numbers by Method")
ax.set_xlabel("Peaks")
ax.set_xlim(0, max(peaks) + 10000)
ax.set_ylabel("Methods")
ax.set_ylim(-0.5, 3.5)

# 调整布局
plt.tight_layout()
plt.show()





