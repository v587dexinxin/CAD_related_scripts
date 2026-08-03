# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 16:34:05 2026

@author: lenovo
"""

import pandas as pd


def read_tad(file):
    """
    Read TAD file with columns: chr, start, end
    """
    df = pd.read_csv(
        file,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        usecols=(0 , 1 , 2)
    )
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def is_same_tad(tad1, tad2, tolerance=50000):
    """
    判断两个 TAD 是否相同：
    chr 相同，且 start 和 end 分别相差不超过 tolerance
    """
    return (
        tad1["chr"] == tad2["chr"]
        and abs(tad1["start"] - tad2["start"]) <= tolerance
        and abs(tad1["end"] - tad2["end"]) <= tolerance
    )


def merge_tad_union(tad_dfs, tolerance=50000):
    """
    对多组 TAD 取并集。
    如果两个 TAD 的左右 boundary 均在 tolerance 范围内，
    则认为是同一个 TAD，只保留一个代表 TAD。
    """
    all_tads = pd.concat(tad_dfs, ignore_index=True)

    # 排序，方便合并
    all_tads = all_tads.sort_values(["chr", "start", "end"]).reset_index(drop=True)

    union_tads = []

    for _, tad in all_tads.iterrows():
        tad_dict = tad.to_dict()

        found = False

        for i, existing_tad in enumerate(union_tads):
            if is_same_tad(tad_dict, existing_tad, tolerance=tolerance):
                found = True

                # 可选：用平均 boundary 作为合并后的代表 TAD
                # union_tads[i]["start"] = int((existing_tad["start"] + tad_dict["start"]) / 2)
                # union_tads[i]["end"] = int((existing_tad["end"] + tad_dict["end"]) / 2)

                # 记录该 union TAD 来自多少个原始 TAD
                union_tads[i]["count"] += 1
                break

        if not found:
            tad_dict["count"] = 1
            union_tads.append(tad_dict)

    union_df = pd.DataFrame(union_tads)
    union_df = union_df.sort_values(["chr", "start", "end"]).reset_index(drop=True)

    return union_df


# =========================
# 输入三个 TAD 文件
# =========================

tad1 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bed")
tad2 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bed")
tad3 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bed")

# 取并集
union_tads = merge_tad_union(
    [tad1, tad2, tad3],
    tolerance=20000
)

union_tads['start1'] = union_tads['start'] - 250
union_tads['end1'] = union_tads['start'] + 250
union_tads['start2'] = union_tads['end'] - 250
union_tads['end2'] = union_tads['end'] + 250

# 保存结果
union_tads[["chr", "start1", "end1" , "chr" , "start2" , "end2"]].to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\WT_LSS_OSS_union_loops.bed",
    sep="\t",
    header=False,
    index=False
)

# 如果想保留 count 信息
union_tads[["chr", "start1", "end1" , "chr" , "start2" , "end2" , "count"]].to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\WT_LSS_OSS_union_loops_with_count.bed",
    sep="\t",
    header=True,
    index=False
)


