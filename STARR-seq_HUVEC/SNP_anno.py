# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 10:32:52 2026

@author: lenovo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
from matplotlib.gridspec import GridSpec

def seqs_to_index(seqs_set):
    seqs_new = []
    for i in seqs_set:
        j = i.split(',')
        seqs_new.extend(j)
        
    seqs_new = list(set(seqs_new))
    seqs_index = SNPs[SNPs['Seq_name'].isin(seqs_new)].index
    
    return seqs_index





SNPs = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)
SNPs['Coatis'] = 0
SNPs['ATAC'] = 0
SNPs['H3K27ac'] = 0
SNPs['eQTL'] = 0
SNPs['enhancer'] = 0
SNPs['silencer'] = 0
SNPs['DAVs'] = 0
SNPs['act_RE'] = 0



Coatis_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional_new\\peaks\\union_peaks\\union_peaks_sorted_merged.narrowPeak' , header = None , usecols = (0 , 1 , 2))
Coatis_peak.columns = ['chr' , 'start' , 'end']
ATAC_peak = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\ATAC-seq\\HUVEC_ATAC_union_peaks_sorted_merged.narrowPeak' , header = None , usecols = (0 , 1 , 2))
ATAC_peak.columns = ['chr' , 'start' , 'end']
H3K27ac_peak = pd.read_table('H:\\work\\literature_data\\HUVEC\\HUVEC_HeK27ac_hg38_ENCFF077LGZ.bed' , header = None , usecols = (0 , 1 , 2))
H3K27ac_peak.columns = ['chr' , 'start' , 'end']
eqtl1 = pd.read_table('H:\\work\\literature_data\\eQTL\\CAD_related\\Artery_Coronary.v11.eGenes.txt')
eqtl1 = eqtl1[eqtl1['qval'] <= 0.05]
eqtl2 = pd.read_table('H:\\work\\literature_data\\eQTL\\CAD_related\\Artery_Tibial.v11.eGenes.txt')
eqtl2 = eqtl2[eqtl2['qval'] <= 0.05]
eqtl3 = pd.read_table('H:\\work\\literature_data\\eQTL\\CAD_related\\Heart_Atrial_Appendage.v11.eGenes.txt')
eqtl3 = eqtl3[eqtl3['qval'] <= 0.05]
eqtl4 = pd.read_table('H:\\work\\literature_data\\eQTL\\CAD_related\\Heart_Left_Ventricle.v11.eGenes.txt')
eqtl4 = eqtl4[eqtl4['qval'] <= 0.05]
eqtl = pd.concat([eqtl1 , eqtl2 , eqtl3 , eqtl4])

eqtl = eqtl.drop_duplicates(subset=['rs_id_dbSNP157_GRCh38p14'])

enhancer_n = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\HUVEC_nonrisk_enhancer_p0.05_fc+-0_edgR.csv' , header = 0)
enhancer_r = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\HUVEC_risk_enhancer_p0.05_fc+-0_edgR.csv' , header = 0)
silencer_n = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\HUVEC_nonrisk_silencer_p0.05_fc+-0_edgR.csv' , header = 0)
silencer_r = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\HUVEC_risk_silencer_p0.05_fc+-0_edgR.csv' , header = 0)



# enhancer_n = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\HUVEC_nonrisk_enhancer_p0.05_fc+-0.3_edgR.csv' , header = 0)
# enhancer_r = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\HUVEC_risk_enhancer_p0.05_fc+-0.3_edgR.csv' , header = 0)
# silencer_n = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\HUVEC_nonrisk_silencer_p0.05_fc+-0.3_edgR.csv' , header = 0)
# silencer_r = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\HUVEC_risk_silencer_p0.05_fc+-0.3_edgR.csv' , header = 0)






############enhancer_silencer
# regulatory_ele = pd.concat([enhancer_n , enhancer_r , silencer_n , silencer_r])    
# regulatory_ele = regulatory_ele.drop_duplicates(subset=['seq'])
enhancer = pd.concat([enhancer_n , enhancer_r])
enhancer = enhancer.drop_duplicates(subset=['seq'])

enhancer_seq = list(set(enhancer['seq']))
enhancer_index = seqs_to_index(enhancer_seq)


silencer = pd.concat([silencer_n , silencer_r])
silencer = silencer.drop_duplicates(subset=['seq'])

silencer_seq = list(set(silencer['seq']))
silencer_index = seqs_to_index(silencer_seq)

############activate_regulatory_elements
act_RE_seq = pd.concat([enhancer , silencer])
act_RE_seq = act_RE_seq.drop_duplicates(subset=['seq'])

act_seq = list(set(act_RE_seq['seq']))
act_index = seqs_to_index(act_seq)




############DAVs
enhancer_speci = set(enhancer_n['seq']) ^ set(enhancer_r['seq'])
silencer_speci = set(silencer_n['seq']) ^ set(silencer_r['seq'])

DAVs_seqs = list(set(enhancer_speci | silencer_speci))
DAVs_index = seqs_to_index(DAVs_seqs)



#####DAVs_ttest
dav_ttest = pd.read_table('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVs\\HUVEC_nonrisk_risk_cDNA_vs_plasmid_ttest.txt' , header = 0)

dav_t_seq = list(set(dav_ttest[dav_ttest['p_value'] <= 0.05]['Geneid']))

DAVs_index_t = seqs_to_index(dav_t_seq)

DAVs_index_final = DAVs_index & DAVs_index_t





for i in SNPs.index:
    seq = SNPs.iloc[i]['Seq_name']
    g = 'chr' + SNPs.iloc[i]['CHR_ID']
    pos_s = SNPs.iloc[i]['CHR_POS'] - 200
    pos_e = SNPs.iloc[i]['CHR_POS'] + 200
    snp = SNPs.iloc[i]['SNPS']
    tmp_coatis = Coatis_peak[Coatis_peak['chr'] == g]
    tmp_atac = ATAC_peak[ATAC_peak['chr'] == g]
    tmp_h3k27ac = H3K27ac_peak[H3K27ac_peak['chr'] == g]
    tmp_e1 = eqtl1[eqtl1['chr'] == g]
    tmp_e2 = eqtl2[eqtl2['chr'] == g]
    tmp_e3 = eqtl3[eqtl3['chr'] == g]
    tmp_e4 = eqtl4[eqtl4['chr'] == g]    

    mask1 = (tmp_coatis['start'] <= pos_e) & (tmp_coatis['end'] >= pos_s)
    mask2 = (tmp_atac['start'] <= pos_e) & (tmp_atac['end'] >= pos_s)
    mask3 = (tmp_h3k27ac['start'] <= pos_e) & (tmp_h3k27ac['end'] >= pos_s)
    mask4 = tmp_e1['rs_id_dbSNP157_GRCh38p14'] == snp
    mask5 = tmp_e2['rs_id_dbSNP157_GRCh38p14'] == snp
    mask6 = tmp_e3['rs_id_dbSNP157_GRCh38p14'] == snp
    mask7 = tmp_e4['rs_id_dbSNP157_GRCh38p14'] == snp

    
    overlap1 = tmp_coatis[mask1]
    overlap2 = tmp_atac[mask2]
    overlap3 = tmp_h3k27ac[mask3]
    overlap4 = tmp_e1[mask4]
    overlap5 = tmp_e2[mask5]
    overlap6 = tmp_e3[mask6]
    overlap7 = tmp_e4[mask7]

    
    overlap_e = pd.concat([overlap1 , overlap2 , overlap3 , overlap4])
    tmp = {'Coatis' : overlap1 , 'ATAC' : overlap2 , 'H3K27ac': overlap3 , 'eQTL' : overlap_e }
    for j in tmp:
        if len(tmp[j]) > 0:
            SNPs.at[i, j] = 1

SNPs.loc[enhancer_index, 'enhancer'] = 1
SNPs.loc[silencer_index, 'silencer'] = 1
SNPs.loc[DAVs_index_final, 'DAVs'] = 1
SNPs.loc[act_index, 'act_RE'] = 1





###########write_to_files
a = SNPs[SNPs['act_RE'] == 1]
a.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVS\\HUVEC_nonrisk_risk_enhancer_silencer_union.csv' , header = True , index = None)

a['start'] = a['CHR_POS'] - 200
a['end'] = a['CHR_POS'] + 200
a['chr'] = 'chr' + a['CHR_ID'].astype(str)

a[['chr' , 'start' , 'end' , 'Seq_name']].to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\fc_0\\DAVS\\HUVEC_nonrisk_risk_enhancer_silencer_union.bed' , header = None , index = None , sep = '\t')



SNPs.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\first6000_SNPs_anno\\SNPs_anno.csv' , header = True , index = None)

# a = SNPs[(SNPs['Coatis'] == 1) & (SNPs['ATAC'] == 1) & (SNPs['H3K27ac'] == 1) & (SNPs['eQTL'] == 1) & ((SNPs['enhancer'] == 1) | (SNPs['silencer'] == 1))]
# a.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\SNPs\\signal_enriched_SNPs.csv' , header = True , index = None)

# b = SNPs[(SNPs['Coatis'] == 1) & (SNPs['eQTL'] == 1) & ((SNPs['enhancer'] == 1) | (SNPs['silencer'] == 1))]
# b.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\SNPs\\Coatis_signal_enriched_SNPs.csv' , header = True , index = None)





# selected_seqs = pd.read_excel('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\verification_experiments\\Total_examples.xlsx')
# selected_seqs_1 = set(selected_seqs['seq'].str.split('_', expand=True)[4])
# selected_seqs_2 = set(selected_seqs['seq'].str.split('_', expand=True)[6])


# selected_seqs = selected_seqs_1 | selected_seqs_2







###################UpSet plots


df = SNPs[(SNPs['Coatis'] == 1)]

anno_cols = ['Coatis', 'eQTL', 'ATAC', 'H3K27ac', 'act_RE', 'DAVs']




# 确保为 0/1 整数
df[anno_cols] = df[anno_cols].fillna(0).astype(int)

# 去掉所有注释都为 0 的 SNP
df_anno = df[df[anno_cols].sum(axis=1) > 0].copy()

# =========================
# 3. 统计每种组合的 SNP 数量
# =========================
combo_counts = (
    df_anno
    .groupby(anno_cols)
    .size()
    .reset_index(name="count")
)

combo_counts = combo_counts[combo_counts["count"] > 0]
combo_counts = combo_counts.sort_values("count", ascending=False).reset_index(drop=True)

# 只展示前 top_n 个组合
top_n = 30
combo_counts = combo_counts.head(top_n)

# 每一类 annotation 的总 SNP 数
set_sizes = df_anno[anno_cols].sum()
set_sizes_ordered = set_sizes.loc[anno_cols]

# =========================
# 4. 作图参数
# =========================
light_blue = "#bdd7e7"   # 右侧横向柱状图淡蓝色
bar_color = "0.45"       # 上方交集柱状图颜色
dot_color = "0.15"       # 实心点颜色
empty_dot_color = "0.85" # 空点颜色

n_combo = combo_counts.shape[0]
n_sets = len(anno_cols)

x = np.arange(n_combo)
y = np.arange(n_sets)

# =========================
# 5. 创建画布
# =========================
fig = plt.figure(figsize=(12, 6))



gs = GridSpec(
    2, 2,
    width_ratios=[5, 1.3],
    height_ratios=[3, 2],
    wspace=0.05,
    hspace=0.05
)

ax_bar = fig.add_subplot(gs[0, 0])
ax_matrix = fig.add_subplot(gs[1, 0], sharex=ax_bar)
ax_setsize = fig.add_subplot(gs[1, 1], sharey=ax_matrix)


# =========================
# 6. 上方交集柱状图
# =========================
ax_bar.bar(
    x,
    combo_counts["count"],
    color=bar_color,
    width=0.7
)

for i, v in enumerate(combo_counts["count"]):
    ax_bar.text(
        i,
        v + max(combo_counts["count"]) * 0.02,
        str(v),
        ha="center",
        va="bottom",
        fontsize=8
    )

ax_bar.set_ylabel("No. of Intersections", fontsize=11)
ax_bar.set_xticks([])
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)

# =========================
# 7. 下方 annotation 组合点阵
#    修改为条纹背景 + 大圆点风格
# =========================

# 背景条纹颜色
row_bg_color = "#f0f0f0"

# 点颜色
active_dot_color = "black"
inactive_dot_color = "#cfcfcf"

# 点大小
active_dot_size = 120
inactive_dot_size = 120

# 先画每一行的浅灰色背景条纹
for j in range(n_sets):
    if j % 2 == 1:
        ax_matrix.axhspan(
            j - 0.5,
            j + 0.5,
            color=row_bg_color,
            zorder=0
        )

# 再画点和连线
for i in range(n_combo):
    values = combo_counts.loc[i, anno_cols].values

    active_y = []

    for j, val in enumerate(values):
        if val == 1:
            ax_matrix.scatter(
                i,
                j,
                s=active_dot_size,
                color=active_dot_color,
                edgecolor=active_dot_color,
                linewidth=0,
                zorder=3
            )
            active_y.append(j)
        else:
            ax_matrix.scatter(
                i,
                j,
                s=inactive_dot_size,
                color=inactive_dot_color,
                edgecolor=inactive_dot_color,
                linewidth=0,
                zorder=2
            )

    # 连接同一组合中为 1 的点
    if len(active_y) >= 2:
        ax_matrix.plot(
            [i, i],
            [min(active_y), max(active_y)],
            color=active_dot_color,
            linewidth=2,
            zorder=1
        )

# 设置 y 轴标签
ax_matrix.set_yticks(y)
ax_matrix.set_yticklabels(
    anno_cols,
    fontsize=11
)

ax_matrix.set_xlabel(
    "Functional annotation combinations",
    fontsize=12
)

ax_matrix.set_xlim(-0.5, n_combo - 0.5)
ax_matrix.set_ylim(-0.5, n_sets - 0.5)

# 是否倒置 y 轴取决于你想显示的顺序
# 如果 anno_cols = ['Coatis', 'eQTL', 'ATAC', 'H3K27ac', 'DAVs', 'act_RE']
# 并希望图中从上到下显示 act_RE, DAVs, H3K27ac, ATAC, eQTL, Coatis，
# 就保留 invert_yaxis()
ax_matrix.invert_yaxis()

# 去掉边框
ax_matrix.spines["top"].set_visible(False)
ax_matrix.spines["right"].set_visible(False)
ax_matrix.spines["bottom"].set_visible(False)
ax_matrix.spines["left"].set_visible(False)

# 去掉 x 轴刻度
ax_matrix.tick_params(
    axis="x",
    bottom=False,
    labelbottom=False
)

# y 轴保留文字，但去掉刻度线
ax_matrix.tick_params(
    axis="y",
    left=False,
    labelleft=True
)



# =========================
# 8. 右侧 set size 横向柱状图
# =========================
ax_setsize.barh(
    y,
    set_sizes_ordered.values,
    color=light_blue,
    height=0.6
)

for j, v in enumerate(set_sizes_ordered.values):
    ax_setsize.text(
        v + max(set_sizes_ordered.values) * 0.02,
        j,
        str(v),
        va="center",
        ha="left",
        fontsize=9
    )

ax_setsize.set_xlabel(
    "Set size",
    fontsize=12
)

ax_setsize.set_yticks(y)

# 关键修改：隐藏右侧 y 轴文字，但不影响左侧共享 y 轴
ax_setsize.tick_params(
    axis="y",
    left=False,
    labelleft=False
)

ax_setsize.invert_yaxis()

ax_setsize.spines["top"].set_visible(False)
ax_setsize.spines["right"].set_visible(False)
ax_setsize.spines["left"].set_visible(False)







plt.savefig(
    "H:\\work\\Postdoctoral\\GWAS疾病位点检测\\论文投稿\\Figures\\Fig4\\Candidate_SNPs_functional_annotation_UpSet_fixed_order.pdf",
    bbox_inches="tight"
)



