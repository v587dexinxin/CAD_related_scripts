# -*- coding: utf-8 -*-
"""
Created on Wed May 20 12:40:08 2026

@author: lenovo
"""



###linux order

bedtools intersect -a ../union_peaks_sorted_merged.bed -b ../HiRPC_WT_allreps_q0.05_peaks_sorted_merged.bed -c | bedtools intersect -a - -b ../HiRPC_LS_allreps_q0.05_peaks_sorted_merged.bed -c | bedtools intersect -a - -b ../HiRPC_OS_allreps_q0.05_peaks_sorted_merged.bed -c > WT_LS_OS.consensus.with_counts.bed
awk 'BEGIN{OFS="\t"}
{
    wt = ($4 > 0)
    ls = ($5 > 0)
    os = ($6 > 0)

    if (wt && !ls && !os) wt_only++
    else if (!wt && ls && !os) ls_only++
    else if (!wt && !ls && os) os_only++
    else if (wt && ls && !os) wt_ls++
    else if (wt && !ls && os) wt_os++
    else if (!wt && ls && os) ls_os++
    else if (wt && ls && os) wt_ls_os++
}
END{
    print "WT_only", wt_only+0
    print "LS_only", ls_only+0
    print "OS_only", os_only+0
    print "WT_LS", wt_ls+0
    print "WT_OS", wt_os+0
    print "LS_OS", ls_os+0
    print "WT_LS_OS", wt_ls_os+0
}' WT_LS_OS.consensus.with_counts.bed > WT_LS_OS.venn_counts.txt

#############################




import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os


################Hi-Coatis
outdir = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\Venn3"
count_file = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\Venn3\\WT_LS_OS.venn_counts.txt"

counts = {}

with open(count_file) as f:
    for line in f:
        key, value = line.strip().split()
        counts[key] = int(value)

# matplotlib-venn 的 subsets 顺序是：
# (100, 010, 110, 001, 101, 011, 111)
# 即：
# WT only, LS only, WT&LS only, OS only, WT&OS only, LS&OS only, WT&LS&OS

subsets = (
    counts["WT_only"],
    counts["LS_only"],
    counts["WT_LS"],
    counts["OS_only"],
    counts["WT_OS"],
    counts["LS_OS"],
    counts["WT_LS_OS"]
)

plt.figure(figsize=(6, 6))

v = venn3(
    subsets=subsets,
    set_labels=("WT", "LS", "OS")
)

plt.title("Overlap of WT, LSS and OSS Hi-Coatis_peaks", fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(outdir, "WT_LS_OS_Hi-Coatis_peak_venn3.pdf"))
# plt.savefig("WT_LS_OS_peak_venn3.png", dpi=300)

plt.close()


################ATAC
outdir = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\Venn3"
count_file = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\peaks\\union_peaks\\Venn3\\WT_LS_OS_ATAC.venn_counts.txt"

counts = {}

with open(count_file) as f:
    for line in f:
        key, value = line.strip().split()
        counts[key] = int(value)

# matplotlib-venn 的 subsets 顺序是：
# (100, 010, 110, 001, 101, 011, 111)
# 即：
# WT only, LS only, WT&LS only, OS only, WT&OS only, LS&OS only, WT&LS&OS

subsets = (
    counts["WT_only"],
    counts["LS_only"],
    counts["WT_LS"],
    counts["OS_only"],
    counts["WT_OS"],
    counts["LS_OS"],
    counts["WT_LS_OS"]
)

plt.figure(figsize=(6, 6))

v = venn3(
    subsets=subsets,
    set_labels=("WT", "LS", "OS")
)

plt.title("Overlap of WT, LS and OS ATAC peaks", fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(outdir, "WT_LS_OS_ATAC_peak_venn3.pdf"))
# plt.savefig("WT_LS_OS_peak_venn3.png", dpi=300)

plt.close()


