# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:05:13 2026

@author: lenovo
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def plot_starr_rank_from_table(
    infile,
    col_activity="log2_cDNA_over_DNA",  # 你的log2(cDNA/DNA)列名
    col_p="pvalue",                   # 你的p值列名
    p_thr=0.05,
    enh_thr=0.5,
    sil_thr=-0.5,
    sep=None,                         # 自动推断分隔符；如是tsv可写"\t"
    title=None,
    outdir=None,
    out_fig=None,
    out_e=None,
    out_s=None
    
):
    # 读取数据
    df = pd.read_csv(infile, sep=sep , header = 0) if sep is not None else pd.read_csv(infile)
    # 如果自动读错分隔符（常见tsv），你可以改成：pd.read_csv(infile, sep="\t")

    # 取需要列并清理
    x = df[[col_activity, col_p]].copy()
    x['index'] = range(len(x))
    x[col_activity] = -pd.to_numeric(x[col_activity], errors="coerce")
    x[col_p] = pd.to_numeric(x[col_p], errors="coerce")
    x = x.dropna(subset=[col_activity, col_p])

    # 按 activity 降序排序并生成 rank
    x = x.sort_values(col_activity, ascending=False).reset_index(drop=True)
    x["Rank"] = np.arange(1, len(x) + 1)

    # 分类
    enh = (x[col_p] < p_thr) & (x[col_activity] > enh_thr)
    sil = (x[col_p] < p_thr) & (x[col_activity] < sil_thr)
    other = ~(enh | sil)

    n_enh = int(enh.sum())
    n_sil = int(sil.sum())

    # 作图
    plt.figure(figsize=(5, 4), dpi=200)

    # 先画灰色，再画彩色覆盖（更像示例图）
    plt.scatter(x.loc[other, "Rank"], x.loc[other, col_activity],
                s=10, c="lightgray", alpha=0.9, edgecolors="none")

    plt.scatter(x.loc[enh, "Rank"], x.loc[enh, col_activity],
                s=10, c="#ff4d4d", alpha=0.9, edgecolors="none",
                label=f"Enhancer: {n_enh}")

    plt.scatter(x.loc[sil, "Rank"], x.loc[sil, col_activity],
                s=10, c="#4d79ff", alpha=0.9, edgecolors="none",
                label=f"Silencer: {n_sil}")

    # 阈值线（可选，方便读图）
    plt.axhline(enh_thr, lw=1, ls="--", color="k", alpha=0.4)
    plt.axhline(sil_thr, lw=1, ls="--", color="k", alpha=0.4)

    plt.xlabel("Rank")
    plt.ylabel("log2(cDNA/DNA)")
    if title:
        plt.title(title)

    plt.legend(frameon=False, loc="upper left")
    plt.tight_layout()
    plt.show()

    if outdir:
        plt.savefig(os.path.join(outdir , out_fig), bbox_inches="tight")
        print("Saved:", os.path.join(outdir , out_fig))
    if out_e:
        index_e = x[enh]['index']
        df.loc[index_e].to_csv(os.path.join(outdir , out_e + '.csv') , header = True , index = None)
        df.loc[index_e][['seqnames' , 'start' , 'end']].to_csv(os.path.join(outdir , out_e + '.bed') , header = None , index = None , sep = '\t')
    if out_s:
        index_s = x[sil]['index']
        df.loc[index_s].to_csv(os.path.join(outdir , out_s + '.csv') , header = True , index = None)
        df.loc[index_s][['seqnames' , 'start' , 'end']].to_csv(os.path.join(outdir , out_s + '.bed') , header = None , index = None , sep = '\t')
                

    

    # 同时返回统计结果方便你写论文/图注
    return {"N_total": len(x), "N_enhancer": n_enh, "N_silencer": n_sil}





################data

HUVEC_nonrisk = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_nonrisk\\HUVEC_nonrisk_cDNA_VS_plasmid_edgeR.csv' , header = 0)
HUVEC_risk = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\HUVEC_risk\\HUVEC_risk_cDNA_VS_plasmid_edgeR.csv' , header = 0)





stats = plot_starr_rank_from_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_nonrisk\\HUVEC_nonrisk_cDNA_VS_plasmid_edgeR.csv',
                                    col_activity="Fold",
                                    col_p="p.value",
                                    p_thr=0.05,
                                    enh_thr=0,
                                    sil_thr=0,
                                    sep=',' , 
                                    title="STARR-seq nonrisk ranked activity",
                                    outdir="H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_nonrisk",
                                    out_fig="HUVEC_nonrisk_starr_rank_p0.05_fc0.pdf",
                                    out_e="HUVEC_nonrisk_cDNA_VS_plasmid_edgeR_enhancer_p0.05_fc0",
                                    out_s="HUVEC_nonrisk_cDNA_VS_plasmid_edgeR_silencer_p0.05_fc0")






stats = plot_starr_rank_from_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\DiffBind\\HUVEC_risk\\HUVEC_risk_cDNA_VS_plasmid_edgeR.csv',
                                    col_activity="Fold",
                                    col_p="p.value",
                                    p_thr=0.05,
                                    enh_thr=0,
                                    sil_thr=0,
                                    sep=',' , 
                                    title="STARR-seq risk ranked activity",
                                    outdir="H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\HUVEC_risk",
                                    out_fig="HUVEC_risk_starr_rank_p0.05_fc0.pdf",
                                    out_e="HUVEC_risk_cDNA_VS_plasmid_edgeR_enhancer_p0.05_fc0",
                                    out_s="HUVEC_risk_cDNA_VS_plasmid_edgeR_silencer_p0.05_fc0")




