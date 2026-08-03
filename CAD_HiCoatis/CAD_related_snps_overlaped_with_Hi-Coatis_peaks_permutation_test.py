# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:40:51 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import subprocess
import os
from tqdm import tqdm
import matplotlib.pyplot as plt


cad_snp_bed = "/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/CAD_SNPs/Hapmap/CAD_related_SNPs_LD0.8_all_+-200bp.bed"
background_snp_bed = "/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/CAD_SNPs/Hapmap/noCAD_related_SNPs_LD0.8_all_background_+-200bp.bed"
peak_bed = "/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/one-dimension_new/mapping_new/macs2/union_peaks/union_peaks_sorted_merged.bed"

outdir = "/scratch/2026-05-25/bio-shenw/Cardiovascular_disease_STARR-seq/CAD_SNPs/CAD_SNPs_VS_Hi-Coatis_peaks/CAD_SNP_chr_matched_random_enrichment_2"
os.makedirs(outdir, exist_ok=True)

n_random = 1000
observed_overlap_count = 1632

random_seed = 123
np.random.seed(random_seed)


cad = pd.read_csv(cad_snp_bed, sep="\t", header=None, comment="#")
bg = pd.read_csv(background_snp_bed, sep="\t", header=None, comment="#")

cad = cad.iloc[:, :4] if cad.shape[1] >= 4 else cad.iloc[:, :3]
bg = bg.iloc[:, :4] if bg.shape[1] >= 4 else bg.iloc[:, :3]

cad.columns = ["chr", "start", "end"] + list(range(3, cad.shape[1]))
bg.columns = ["chr", "start", "end"] + list(range(3, bg.shape[1]))

n_cad_snps = cad.shape[0]
observed_overlap_ratio = observed_overlap_count / n_cad_snps

cad_chr_counts = cad["chr"].value_counts().to_dict()

random_results = []

for i in tqdm(range(n_random)):

    sampled_list = []

    for chrom, n_chr in cad_chr_counts.items():

        bg_chr = bg[bg["chr"] == chrom]

        if bg_chr.shape[0] < n_chr:
            raise ValueError(
                f"Not enough background SNPs on {chrom}: "
                f"need {n_chr}, available {bg_chr.shape[0]}"
            )

        sampled_chr = bg_chr.sample(
            n=n_chr,
            replace=False,
            random_state=random_seed + i
        )

        sampled_list.append(sampled_chr)

    sampled_bg = pd.concat(sampled_list, axis=0)

    random_bed = f"{outdir}/random_chr_matched_{i+1}.bed"
    overlap_bed = f"{outdir}/random_chr_matched_{i+1}_overlap.bed"

    sampled_bg.to_csv(random_bed, sep="\t", header=False, index=False)

    cmd = [
        "bedtools", "intersect",
        "-a", random_bed,
        "-b", peak_bed,
        "-wa", "-u"
    ]

    with open(overlap_bed, "w") as fout:
        subprocess.run(cmd, stdout=fout, check=True)

    overlap_count = sum(1 for _ in open(overlap_bed))
    overlap_ratio = overlap_count / n_cad_snps

    random_results.append({
        "iteration": i + 1,
        "overlap_count": overlap_count,
        "overlap_ratio": overlap_ratio
    })

    # os.remove(random_bed)
    # os.remove(overlap_bed)


random_df = pd.DataFrame(random_results)
random_df.to_csv(
    f"{outdir}/chr_matched_random_overlap_results.tsv",
    sep="\t",
    index=False
)

random_mean = random_df["overlap_ratio"].mean()
random_sd = random_df["overlap_ratio"].std()
fold_enrichment = observed_overlap_ratio / random_mean

empirical_p = (
    (random_df["overlap_ratio"] >= observed_overlap_ratio).sum() + 1
) / (n_random + 1)

summary = pd.DataFrame({
    "observed_overlap_count": [observed_overlap_count],
    "observed_total_snps": [n_cad_snps],
    "observed_overlap_ratio": [observed_overlap_ratio],
    "random_mean_overlap_ratio": [random_mean],
    "random_sd_overlap_ratio": [random_sd],
    "fold_enrichment": [fold_enrichment],
    "empirical_p_value": [empirical_p],
    "n_random": [n_random]
})

summary.to_csv(
    f"{outdir}/chr_matched_random_enrichment_summary.tsv",
    sep="\t",
    index=False
)

print(summary.T)



########


random_df = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\CAD_SNPs_VS_Hi-Coatis_peaks\\CAD_SNP_chr_matched_random_enrichment_2\\chr_matched_random_overlap_results.tsv' , header = 0)

result = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\CAD_SNPs_VS_Hi-Coatis_peaks\\CAD_SNP_chr_matched_random_enrichment_2\\chr_matched_random_enrichment_summary.tsv' , header = 0)


observed_overlap_ratio = result['observed_overlap_ratio'].iloc[0]
fold_enrichment = result['fold_enrichment'].iloc[0]
empirical_p = result['empirical_p_value'].iloc[0]



plt.figure(figsize=(5.5, 4.5))

plt.hist(
    random_df["overlap_ratio"] * 100,
    bins=30,
    alpha=0.75,
    edgecolor="black"
)

plt.axvline(
    observed_overlap_ratio * 100,
    color="red",
    linestyle="--",
    linewidth=2,
    label=f"CAD SNPs = {observed_overlap_ratio * 100:.2f}%"
)

plt.xlabel("Overlap ratio with HUVEC Hi-Coatis peaks (%)")
plt.ylabel("Random sampling frequency")
plt.title("Chromosome-matched random background")

plt.legend(frameon=False)

plt.text(
    0.98,
    0.95,
    f"Fold enrichment = {fold_enrichment:.2f}\nEmpirical P = {empirical_p:.3g}",
    transform=plt.gca().transAxes,
    ha="right",
    va="top"
)

plt.tight_layout()
plt.savefig(
    f"{outdir}/CAD_SNP_chr_matched_random_enrichment_1.pdf",
    dpi=300
)
plt.show()