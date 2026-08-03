# -*- coding: utf-8 -*-
"""
Created on Wed May 13 18:09:14 2026

@author: lenovo
"""

import pandas as pd





def classify_peak(annotation):
    annotation = str(annotation)

    if "Promoter" in annotation:
        return "Promoter-associated"

    elif "Distal Intergenic" in annotation:
        return "Distal enhancer-like candidate"

    elif "Intron" in annotation:
        return "Intronic enhancer-like candidate"

    elif "Exon" in annotation:
        return "Exonic regulatory region"

    elif "UTR" in annotation:
        return "UTR-associated regulatory region"

    else:
        return "Other"


def is_candidate_enhancer(feature_class):
    if feature_class in [
        "Distal enhancer-like candidate",
        "Intronic enhancer-like candidate"
    ]:
        return True
    else:
        return False


def crispri_priority(feature_class):
    if feature_class == "Distal enhancer-like candidate":
        return "High"

    elif feature_class == "Intronic enhancer-like candidate":
        return "High"

    elif feature_class in [
        "UTR-associated regulatory region",
        "Exonic regulatory region"
    ]:
        return "Medium"

    else:
        return "Low"





os_speci_peaks = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\specific_peaks_overlaped_ATAC_Anno\\HUVEC_OS_specific_peaks_edgR_q0.05_fc1_overlaped_ATAC_Anno.csv' , header = 0)


df = os_speci_peaks.copy()

# =========================
# 3. 添加新分类列
# =========================
df["feature_class"] = df["annotation"].apply(classify_peak)
df["is_candidate_enhancer"] = df["feature_class"].apply(is_candidate_enhancer)
df["crispri_priority"] = df["feature_class"].apply(crispri_priority)

# =========================
# 4. 输出分类后的文件
# =========================
df[['seqnames', 'start', 'end', 'width', 'strand', 'annotation','feature_class', 'is_candidate_enhancer', 'crispri_priority']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\DiffBind\\specific_peaks_edgR\\specific_peaks_overlaped_ATAC_Anno\\HUVEC_OS_specific_peaks_edgR_q0.05_fc1_overlaped_ATAC_Anno_classify.csv', index=False)

print("Classification finished!")
print(df["feature_class"].value_counts())
print("\nCandidate enhancer number:")
print(df["is_candidate_enhancer"].value_counts())



df[df['']]



