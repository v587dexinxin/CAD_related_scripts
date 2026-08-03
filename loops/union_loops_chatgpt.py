# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 15:50:26 2026

@author: lenovo
"""

import pandas as pd
from collections import defaultdict

# =========================
# 参数设置
# =========================
WT_FILE = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_WT_HiCoatis_oneanchor_binding_loops.bedpe"
LSS_FILE = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_LS_HiCoatis_oneanchor_binding_loops.bedpe"
OSS_FILE = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\HUVEC_OS_HiCoatis_oneanchor_binding_loops.bedpe"

OUTPUT_FILE = "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\loops\\WT_LSS_OSS_union_loops_chatgpt.bedpe"

DISTANCE = 20_000
BIN_SIZE = DISTANCE

columns = [
    "chr1", "start1", "end1",
    "chr2", "start2", "end2"
]


def read_loop_file(filename, condition):
    """读取 BEDPE loop 文件。"""
    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        comment="#",
        usecols=range(6),
        names=columns
    )

    for col in ["start1", "end1", "start2", "end2"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna().copy()

    for col in ["start1", "end1", "start2", "end2"]:
        df[col] = df[col].astype(int)

    df["condition"] = condition
    return df


def normalize_anchor_order(row):
    """
    确保 anchor1 排在 anchor2 前面。

    对跨染色体 loop，按照染色体名称和坐标排序；
    对同染色体 loop，按照中心位置排序。
    """
    center1 = (row["start1"] + row["end1"]) // 2
    center2 = (row["start2"] + row["end2"]) // 2

    anchor1 = (
        row["chr1"],
        row["start1"],
        row["end1"],
        center1
    )

    anchor2 = (
        row["chr2"],
        row["start2"],
        row["end2"],
        center2
    )

    if (anchor1[0], anchor1[3]) > (anchor2[0], anchor2[3]):
        row["chr1"], row["chr2"] = row["chr2"], row["chr1"]
        row["start1"], row["start2"] = row["start2"], row["start1"]
        row["end1"], row["end2"] = row["end2"], row["end1"]

    return row


# =========================
# 读取文件
# WT最先，因此重复时优先保留WT；
# 随后依次保留LSS和OSS中首次出现的loop
# =========================
wt = read_loop_file(WT_FILE, "WT")
lss = read_loop_file(LSS_FILE, "LSS")
oss = read_loop_file(OSS_FILE, "OSS")

all_loops = pd.concat(
    [wt, lss, oss],
    axis=0,
    ignore_index=True
)

all_loops = all_loops.apply(normalize_anchor_order, axis=1)

# 保存最终保留的loop
union_rows = []

# 空间索引：
# key = (chr1, chr2, anchor1_bin, anchor2_bin)
# value = 已保留loop在union_rows中的索引
loop_index = defaultdict(list)


for _, row in all_loops.iterrows():

    chr1 = row["chr1"]
    chr2 = row["chr2"]

    center1 = (row["start1"] + row["end1"]) // 2
    center2 = (row["start2"] + row["end2"]) // 2

    bin1 = center1 // BIN_SIZE
    bin2 = center2 // BIN_SIZE

    is_duplicate = False

    # 距离小于20 kb的loop只可能位于当前bin或相邻bin
    for offset1 in [-1, 0, 1]:
        for offset2 in [-1, 0, 1]:

            key = (
                chr1,
                chr2,
                bin1 + offset1,
                bin2 + offset2
            )

            for retained_index in loop_index.get(key, []):
                retained = union_rows[retained_index]

                retained_center1 = (
                    retained["start1"] + retained["end1"]
                ) // 2

                retained_center2 = (
                    retained["start2"] + retained["end2"]
                ) // 2

                left_distance = abs(center1 - retained_center1)
                right_distance = abs(center2 - retained_center2)

                # 左右两个anchor均满足距离条件
                if (
                    left_distance < DISTANCE
                    and right_distance < DISTANCE
                ):
                    is_duplicate = True
                    break

            if is_duplicate:
                break

        if is_duplicate:
            break

    # 没有匹配到已有loop，则加入并集
    if not is_duplicate:
        retained_index = len(union_rows)

        union_rows.append({
            "chr1": chr1,
            "start1": int(row["start1"]),
            "end1": int(row["end1"]),
            "chr2": chr2,
            "start2": int(row["start2"]),
            "end2": int(row["end2"]),
            "condition": row["condition"]
        })

        key = (
            chr1,
            chr2,
            bin1,
            bin2
        )

        loop_index[key].append(retained_index)


# =========================
# 输出结果
# =========================
union_df = pd.DataFrame(union_rows)

union_df[columns].to_csv(
    OUTPUT_FILE,
    sep="\t",
    header=False,
    index=False
)

print("WT loop数量：", len(wt))
print("LSS loop数量：", len(lss))
print("OSS loop数量：", len(oss))
print("合并前loop总数：", len(all_loops))
print("并集loop数量：", len(union_df))
print("结果文件：", OUTPUT_FILE)

print("\n并集中各条件首次贡献的loop数量：")
print(union_df["condition"].value_counts())