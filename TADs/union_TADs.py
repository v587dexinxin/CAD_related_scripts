# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 20:24:26 2026

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

tad1 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\HUVEC_WT_merged_hg38_balanced_25K_domains.bed")
tad2 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\HUVEC_LS_merged_hg38_balanced_25K_domains.bed")
tad3 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\HUVEC_OS_merged_hg38_balanced_25K_domains.bed")

# 取并集
union_tads = merge_tad_union(
    [tad1, tad2, tad3],
    tolerance=50000
)

# 保存结果
union_tads[["chr", "start", "end"]].to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\union_TADs_25K_domains.bed",
    sep="\t",
    header=False,
    index=False
)

# 如果想保留 count 信息
union_tads.to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\union_TADs_25K_domains_with_count.bed",
    sep="\t",
    header=True,
    index=False
)


##############TAD to boundary
# 左右边界合并
left = union_tads[["chr", "start"]].rename(columns={"start": "boundary"})
right = union_tads[["chr", "end"]].rename(columns={"end": "boundary"})

union_boundary = pd.concat([left, right], ignore_index=True)

# 去重并排序
union_boundary = (
    union_boundary
    .drop_duplicates()
    .sort_values(["chr", "boundary"])
    .reset_index(drop=True)
)

# 同一染色体上，距离小于 50 kb 时保留前面的 boundary
filtered_rows = []

for chrom, group in union_boundary.groupby("chr", sort=False):
    group = group.sort_values("boundary")

    last_kept = None

    for _, row in group.iterrows():
        current = row["boundary"]

        if last_kept is None or current - last_kept >= 50000:
            filtered_rows.append(row)
            last_kept = current

union_boundary_filtered = pd.DataFrame(filtered_rows).reset_index(drop=True)

union_boundary_filtered.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_25K\\union_boundary_25K.bed' , header = None , index = None , sep = '\t')







# =========================
# 输入三个 TAD 文件_10K
# =========================

tad1 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\HUVEC_WT_merged_hg38_balanced_10K_domains.bed")
tad2 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\HUVEC_LS_merged_hg38_balanced_10K_domains.bed")
tad3 = read_tad("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\HUVEC_OS_merged_hg38_balanced_10K_domains.bed")

# 取并集
union_tads = merge_tad_union(
    [tad1, tad2, tad3],
    tolerance=50000
)

# 保存结果
union_tads[["chr", "start", "end"]].to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\union_TADs_10K_domains.bed",
    sep="\t",
    header=False,
    index=False
)


union_tads[["chr", "start", "end" , "chr", "start", "end"]].to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\union_TADs_10K_domains.bedpe",
    sep="\t",
    header=False,
    index=False
)



# 如果想保留 count 信息
union_tads.to_csv(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\union_TADs_10K_domains_with_count.bed",
    sep="\t",
    header=True,
    index=False
)



##############TAD to boundary
# 左右边界合并
left = union_tads[["chr", "start"]].rename(columns={"start": "boundary"})
right = union_tads[["chr", "end"]].rename(columns={"end": "boundary"})

union_boundary = pd.concat([left, right], ignore_index=True)

# 去重并排序
union_boundary = (
    union_boundary
    .drop_duplicates()
    .sort_values(["chr", "boundary"])
    .reset_index(drop=True)
)

# 同一染色体上，距离小于 50 kb 时保留前面的 boundary
filtered_rows = []

for chrom, group in union_boundary.groupby("chr", sort=False):
    group = group.sort_values("boundary")

    last_kept = None

    for _, row in group.iterrows():
        current = row["boundary"]

        if last_kept is None or current - last_kept >= 50000:
            filtered_rows.append(row)
            last_kept = current

union_boundary_filtered = pd.DataFrame(filtered_rows).reset_index(drop=True)

union_boundary_filtered.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\TADs\\res_10K\\union_boundary_10K.bed' , header = None , index = None , sep = '\t')




