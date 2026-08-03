# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:06:59 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    
def pos_overlap_peaks_6000(pos , peaks):

    n = 0 ; data_new = []
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    for g in chrom:
        # print (g)
        tmp_data = pos[pos['CHR_ID'] == g.lstrip('chr')]
        tmp_peaks = peaks[peaks['chr'] == g]
        for i in tmp_data.index:
            p = int(tmp_data.loc[i]['CHR_POS'])
            start = p - 200
            end = p + 200
            seq = tmp_data.loc[i]['Seq_name']
            rs = tmp_data.loc[i]['SNPS']
            mask = (tmp_peaks['start'] <= end) & (tmp_peaks['end'] >= start)
            overlap = tmp_peaks[mask]
            if len(overlap) > 0:
                n += 1
                data_new.append((g , p , rs , seq))
    data_new = pd.DataFrame(data_new , columns=['chr' , 'pos' , 'SNPs' , 'seq_name'])
    data_new = data_new.drop_duplicates(subset=['seq_name'] , keep = 'first')
    print (n)
    return(data_new)
        


def pos_overlap_peaks_1(pos , peaks):

    n = 0
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    peaks_new = []
    for g in chrom:
        # print (g)
        tmp_data = pos[pos['chr'] == g]
        tmp_peaks = peaks[peaks['chr'] == g]
        for i in tmp_data.index:
            start = int(tmp_data.loc[i]['pos']) - 200
            end = int(tmp_data.loc[i]['pos']) + 200
            mask = (tmp_peaks['start'] <= end) & (tmp_peaks['end'] >= start)
            overlap = tmp_peaks[mask]
            if len(overlap) > 0:
                peaks_new.append((g , start , end))
                n += 1
                # print (tmp_data.loc[i]['chr'])
    peaks_new = pd.DataFrame(peaks_new)
    peaks_new = peaks_new.drop_duplicates(peaks_new)
    
    print (n)
    return(peaks_new)
        









#############peaks
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
peaks_a = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\union_peaks_sorted_merged.bed' , header = None)
peaks_a.columns = ['chr' , 'start' , 'end']
peaks_a = peaks_a[peaks_a['chr'].isin(chrom)]
peaks_a.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\union_peaks_sorted_merged_clean.bed' , header = None , index = None , sep = '\t')

# peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\LS_OS_union_peaks_sorted_merged.bed' , header = None)
# peaks.columns = ['chr' , 'start' , 'end']

peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_common_peaks_q0.05_fc0.5_clean.bed' , header = None)
peaks.columns = ['chr' , 'start' , 'end']

# peaks = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\HiRPC_OS_allreps_q0.05_peaks_sorted_merged.bed' , header = None)
# peaks.columns = ['chr' , 'start' , 'end']


# peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.bed' , header = None)
# peaks.columns = ['chr' , 'start' , 'end']


# peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.bed' , header = None)
# peaks.columns = ['chr' , 'start' , 'end']



# deseq2_pos = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_deseq2_RLE_BACKGROUND_minOverlap1_all.csv' , header=0)

# deseq2_pos.columns = ['peak_id', 'Chr', 'Start', 'End', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'direction']

# diff = deseq2_pos[deseq2_pos['padj'] <= 0.05]

# oss_speci = diff[diff['log2FoldChange'] >= 0.5]
# lss_speci = diff[diff['log2FoldChange'] < -0.5]

# lss_speci.columns = ['peak_id', 'chr', 'start', 'end', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'direction']
# oss_speci.columns = ['peak_id', 'chr', 'start', 'end', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'direction']





lss_speci =  pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_LSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
lss_speci.columns = ['chr' , 'start' , 'end']
oss_speci =  pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_OSS_specific_peaks_q0.05_fc0.5_clean.bed' , header = None)
oss_speci.columns = ['chr' , 'start' , 'end']



# wt =  pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\Diffbind_DEseq2\\HUVEC_LS_VS_OS_DESeq2_common_peaks_q0.05_fc0.5_clean.bed' , header = None)





#########CAD_related_pos_first6000
CAD = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)




snps_a = pos_overlap_peaks_6000(CAD , peaks_a)



snps_a.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_SNPs_VS_Hi-Coatis\\CAD_SNPs_first6000_overlaped_with_union_Hi-Coatis_peaks.csv' , header = True , index = None)






snps_c = pos_overlap_peaks_6000(CAD , peaks)
snps_ls = pos_overlap_peaks_6000(CAD , lss_speci)
snps_os = pos_overlap_peaks_6000(CAD , oss_speci)





############CAD_related_Snps
CAD_realated_snps = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\Hapmap\\CAD_related_SNPs_LD0.8_all.bed' , header = None)

CAD_realated_snps.columns = ['chr' , 'pos' , 'SNPs']






overlap_c = pos_overlap_peaks_1(CAD_realated_snps , peaks)
overlap_ls = pos_overlap_peaks_1(CAD_realated_snps , lss_speci)
overlap_os = pos_overlap_peaks_1(CAD_realated_snps , oss_speci)
overlap_u = pos_overlap_peaks_1(CAD_realated_snps , peaks_a)


overlap_c.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\CAD_SNPs_VS_LS_OS_common_Hi-Coatis_peaks\\CAD_SNPs_LD0.8_VS_LS_OS_common_Hi-Coatis_peaks_snp_+-200bp.bed' , header = None , index = None , sep = '\t')
overlap_ls.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\CAD_SNPs_VS_LSS_speci_Hi-Coatis_peaks\\CAD_SNPs_LD0.8_VS_LSS_speci_Hi-Coatis_peaks_snp_+-200bp.bed' , header = None , index = None , sep = '\t')
overlap_os.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\CAD_SNPs_VS_OSS_speci_Hi-Coatis_peaks\\CAD_SNPs_LD0.8_VS_OSS_speci_Hi-Coatis_peaks_snp_+-200bp.bed' , header = None , index = None , sep = '\t')
overlap_u.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\CAD_SNPs_VS_union_Hi-Coatis_peaks\\CAD_SNPs_LD0.8_VS_union_Hi-Coatis_peaks_snp_+-200bp.bed' , header = None , index = None , sep = '\t')



import pandas as pd
import matplotlib.pyplot as plt

# =========================
# 1. 输入数据
# =========================

data = {
    "Condition": ["common", "LSS_speci", "OSS_speci"],
    "Overlap_ratio": [len(overlap_c) / len(peaks), len(overlap_ls) / len(lss_speci), len(overlap_os) / len(oss_speci)]
}

df = pd.DataFrame(data)

# 转换为百分比
df["Overlap_percentage"] = df["Overlap_ratio"] * 100

print(df)

# =========================
# 2. 画柱状图
# =========================

plt.figure(figsize=(4.2, 4.5))

colors = {
    "common": "#4C78A8",    # 蓝色
    "LSS_speci": "#59A14F",   # 绿色
    "OSS_speci": "#D37295"    # 淡紫红色
}

bars = plt.bar(
    df["Condition"],
    df["Overlap_percentage"],
    color=[colors[x] for x in df["Condition"]],
    width=0.6,
    edgecolor="black",
    linewidth=0.8
)

# =========================
# 3. 添加数值标签
# =========================

for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        height + 0.05,
        f"{height:.2f}%",
        ha="center",
        va="bottom",
        fontsize=11
    )

# =========================
# 4. 美化图形
# =========================

plt.ylabel("CAD SNP-overlapped peaks (%)", fontsize=12)
plt.xlabel("")
plt.title("Overlap of CAD SNPs with Hi-Coatis peaks", fontsize=13)

plt.ylim(0, max(df["Overlap_percentage"]) * 1.35)

plt.xticks(fontsize=12)
plt.yticks(fontsize=11)

# 去掉上边框和右边框
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

# =========================
# 5. 保存图片
# =========================

plt.savefig("CAD_SNP_overlap_HiCoatis_peaks_barplot.pdf", dpi=300, bbox_inches="tight")
plt.savefig("CAD_SNP_overlap_HiCoatis_peaks_barplot.png", dpi=300, bbox_inches="tight")

plt.show()
















