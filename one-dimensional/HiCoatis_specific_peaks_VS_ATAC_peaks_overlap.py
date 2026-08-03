# -*- coding: utf-8 -*-
"""
Created on Wed May 13 17:14:14 2026

@author: lenovo
"""



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def Common_peaks(peaks1 , peaks2):
    n = 0 ; common = []
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    
    for g in chrom:
        tmp1 = peaks1[peaks1['chr'] == g]
        tmp2 = peaks2[peaks2['chr'] == g]
        for i in tmp1.index:
            start = tmp1.loc[i]['start']
            end = tmp1.loc[i]['end']
            mask = (tmp2['start'] <= end) & (tmp2['end'] >= start)
            overlap = tmp2[mask]
            if len(overlap) != 0:
                n += 1
                c = tuple(pd.concat([tmp1.loc[i] , overlap.iloc[0]] , axis = 0))
                common.append(c)
    print (n)
    common = pd.DataFrame(common)
    
    return common




##############specific Hi-Coatis_VS_ATAC
speci_ls = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
speci_ls.columns = ['chr' , 'start' , 'end']
speci_os = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
speci_os.columns = ['chr' , 'start' , 'end']



ls_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_LS_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
ls_atac_peak.columns = ['chr' , 'start' , 'end']
os_atac_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\HUVEC_ATAC_OS_peaks_sorted_merged.narrowPeak'  , header = None , usecols=(0 , 1 , 2))
os_atac_peak.columns = ['chr' , 'start' , 'end']



speci_ls_common = Common_peaks(speci_ls , ls_atac_peak)
speci_os_common = Common_peaks(speci_os , os_atac_peak)

speci_ls_common[[0 , 1 , 2]].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\overlaped_with_ATAC\\LSS_DiffBind_DEseq2_specific_peaks_overlaped_ATAC_peaks.bed' , header = None , index = None , sep = '\t')
speci_os_common[[0 , 1 , 2]].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\overlaped_with_ATAC\\OSS_DiffBind_DEseq2_specific_peaks_overlaped_ATAC_peaks.bed' , header = None , index = None , sep = '\t')




# =========================
# 1. 输入数据
# =========================
data = pd.DataFrame({
    "Group": ["LSS-specific\nHi-Coatis peaks", "OSS-specific\nHi-Coatis peaks"],
    "Overlap_with_ATAC": [len(speci_ls_common), len(speci_os_common)],
    "Not_overlap_with_ATAC": [len(speci_ls) - len(speci_ls_common), len(speci_os) - len(speci_os_common)]
})

# 计算总数和比例
data["Total"] = data["Overlap_with_ATAC"] + data["Not_overlap_with_ATAC"]
data["Overlap_ratio"] = data["Overlap_with_ATAC"] / data["Total"] * 100
data["Not_overlap_ratio"] = data["Not_overlap_with_ATAC"] / data["Total"] * 100

print(data)

# =========================
# 2. 画比例堆叠柱状图
# =========================
x = np.arange(len(data))

fig, ax = plt.subplots(figsize=(5, 5))

# 未重叠部分
bar1 = ax.bar(
    x,
    data["Not_overlap_ratio"],
    label="Not overlapped with ATAC",
    color="#c7c7c7",
    edgecolor="black",
    linewidth=0.8
)

# 重叠部分
bar2 = ax.bar(
    x,
    data["Overlap_ratio"],
    bottom=data["Not_overlap_ratio"],
    label="Overlapped with ATAC",
    color="#87CEFA",
    edgecolor="black",
    linewidth=0.8
)

# =========================
# 3. 添加百分比和 peak 数标注
# =========================
for i in range(len(data)):
    # not overlap 标注
    ax.text(
        x[i],
        data.loc[i, "Not_overlap_ratio"] / 2,
        f'{data.loc[i, "Not_overlap_with_ATAC"]}\n'
        f'({data.loc[i, "Not_overlap_ratio"]:.1f}%)',
        ha="center",
        va="center",
        fontsize=10
    )
    
    # overlap 标注
    ax.text(
        x[i],
        data.loc[i, "Not_overlap_ratio"] + data.loc[i, "Overlap_ratio"] / 2,
        f'{data.loc[i, "Overlap_with_ATAC"]}\n'
        f'({data.loc[i, "Overlap_ratio"]:.1f}%)',
        ha="center",
        va="center",
        fontsize=10
    )

# =========================
# 4. 美化图形
# =========================
ax.set_xticks(x)
ax.set_xticklabels(data["Group"], fontsize=11)
ax.set_ylabel("Percentage of Hi-Coatis-specific peaks (%)", fontsize=12)
ax.set_ylim(0, 100)

ax.set_title(
    "Overlap of Hi-Coatis-specific peaks with ATAC-seq peaks",
    fontsize=13
)

ax.legend(
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.15),
    ncol=1,
    fontsize=10
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()


# 保存
plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig3\\HiCoatis_specific_peaks_ATAC_overlap_stacked_bar.pdf")


plt.show()


