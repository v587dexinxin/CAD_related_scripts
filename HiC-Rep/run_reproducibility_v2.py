# -*- coding: utf-8 -*-
"""
Created on Wed May  6 19:18:44 2026

@author: lenovo
"""

from __future__ import division
import sys
import scipy
import numpy as np
import pandas as pd
import cooler
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from palettable.colorbrewer.qualitative import  Paired_10

# colors = Paired_10.hex_colors

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')


def Load_matrix(cooler_fil):
    
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    c = cooler.Cooler(cooler_fil )
    Matrix = {}
    
    for g in chrom:
        Matrix[g] = c.matrix(balance=True).fetch(g)

    
    return Matrix


def  get_Laplacian(M):
     S=M.sum(1)
     i_nz=np.where(S>0)[0]
     S=S[i_nz]
     M=(M[i_nz].T)[i_nz].T
     S=1/np.sqrt(S)
     M=S*M
     M=(S*M.T).T
     n=np.size(S)
     M=np.identity(n)-M
     M=(M+M.T)/2
     return M

def evec_distance(v1,v2):
    d1=np.dot(v1-v2,v1-v2)
    d2=np.dot(v1+v2,v1+v2)
    if d1<d2:
         d=d1
    else:
        d=d2
    return np.sqrt(d)

def get_ipr(evec):
      ipr=1.0/(evec*evec*evec*evec).sum()
      return ipr


def get_reproducibility(M1,M2,num_evec):
   k1=np.sign(M1).sum(1)
   d1=np.diag(M1)
   kd1=~((k1==1)*(d1>0))
   k2=np.sign(M2).sum(1)
   d2=np.diag(M2)
   kd2=~((k2==1)*(d2>0))
   iz=np.nonzero((k1+k2>0)*(kd1>0)*(kd2>0))[0]
   M1b=(M1[iz].T)[iz].T
   M2b=(M2[iz].T)[iz].T

   i_nz1=np.where(M1b.sum(1)>0)[0]
   i_nz2=np.where(M2b.sum(1)>0)[0]

   
   M1b_L=get_Laplacian(M1b)
   M2b_L=get_Laplacian(M2b)
   
   a1, b1=eigsh(M1b_L,k=num_evec,which="SM")
   a2, b2=eigsh(M2b_L,k=num_evec,which="SM")
   
   b1_extend=np.zeros((np.size(M1b,0),num_evec))
   b2_extend=np.zeros((np.size(M2b,0),num_evec))
   for i in range(num_evec):
       b1_extend[i_nz1,i]=b1[:,i]
       b2_extend[i_nz2,i]=b2[:,i]
   
   ipr_cut=5
   ipr1=np.zeros(num_evec)
   ipr2=np.zeros(num_evec)
   for i in range(num_evec):
       ipr1[i]=get_ipr(b1_extend[:,i])
       ipr2[i]=get_ipr(b2_extend[:,i])
  
   b1_extend_eff=b1_extend[:,ipr1>ipr_cut]
   b2_extend_eff=b2_extend[:,ipr2>ipr_cut]
   num_evec_eff=min(np.size(b1_extend_eff,1),np.size(b2_extend_eff,1))
  
   evd=np.zeros(num_evec_eff)
   for i in range(num_evec_eff):
       evd[i]=evec_distance(b1_extend_eff[:,i],b2_extend_eff[:,i])
   
   Sd=evd.sum()
   l=np.sqrt(2)
   evs=abs(l-Sd/num_evec_eff)/l

   N=float(M1.shape[1]);
   if (np.sum(ipr1>N/100)<=1)|(np.sum(ipr2>N/100)<=1):
      print("at least one of the maps does not look like typical Hi-C maps")
      evs = 0
      return evs
   else:
      print("size of maps: %d" %(np.size(M1,0)))
      print("reproducibility score: %6.3f " %(evs))
      print("num_evec_eff: %d" %(num_evec_eff))
      return evs





def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    

#--------------------------------------K562_HIRPC_RNAPII_insituHiC--------------------



matrix_wt = Load_matrix('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_WT_merged.hg38.nodups.mapq_30.1000.mcool::/resolutions/500000' )
matrix_ls = Load_matrix('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_LS_merged.hg38.nodups.mapq_30.1000.mcool::/resolutions/500000' )
matrix_os = Load_matrix('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/merged_all_reps/cool/HUVEC_OS_merged.hg38.nodups.mapq_30.1000.mcool::/resolutions/500000')



Matrix = {'WT' : matrix_wt , 'LSS' : matrix_ls , 'OSS' : matrix_os}



chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
num_evec=5
all_reps = {}



for c in ['WT_LSS' , 'WT_OSS' , 'LSS_OSS']:
    R1 = c.split('_')[0]
    R2 = c.split('_')[1]
    x = []
    
    for g in chrom:
        matrix1 = Matrix[R1][g]
        matrix2 = Matrix[R2][g]
        matrix1[np.isnan(matrix1)] = 0
        matrix2[np.isnan(matrix2)] = 0        
        r = get_reproducibility(matrix1 , matrix2 , num_evec)
        x.append(float(r))
    all_reps[c] = x


############old style
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)

ax.boxplot(all_reps['WT_LSS'] , positions=[1] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'SeaGreen','linewidth':2},
        medianprops={'color':'SeaGreen','linewidth':2},
        capprops={'color':'SeaGreen','linewidth':2},
        whiskerprops={'color':'SeaGreen','linewidth':2})
ax.boxplot(all_reps['WT_OSS'] , positions=[2] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'Goldenrod','linewidth':2},
        medianprops={'color':'Goldenrod','linewidth':2},
        capprops={'color':'Goldenrod','linewidth':2},
        whiskerprops={'color':'Goldenrod','linewidth':2})
ax.boxplot(all_reps['LSS_OSS'] , positions=[3] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'MediumPurple','linewidth':2},
        medianprops={'color':'MediumPurple','linewidth':2},
        capprops={'color':'MediumPurple','linewidth':2},
        whiskerprops={'color':'MediumPurple','linewidth':2})


ax.set_xticks([1 , 2 , 3])
ax.set_xticklabels(['WT_VS_LSS' , 'WT_VS_OSS' , 'LSS_VS_OSS'] , fontsize = 28 , rotation=45)
ax.set_xlim((0.5 , 3.5))
ax.set_ylim((0.2 , 1))
    
ax.set_ylabel('Correlation' , fontsize = 30)
plt.tight_layout()

plt.title('Correlation between WT_VS_LSS_VS_OSS' , fontsize = 20)



run_Plot(fig , '/scratch/2026-05-02/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/HiC_Rep/HUVEC_WT_LSS_OSS_reproducibility_boxplot.pdf')







#########new_style
comparisons = ['WT_LSS' , 'WT_OSS' , 'LSS_OSS']
plot_data = [all_reps[k] for k in ['WT_LSS' , 'WT_OSS' , 'LSS_OSS']]


fig, ax = plt.subplots(figsize=(6, 5))

bp = ax.boxplot(
    plot_data,
    labels=comparisons,
    patch_artist=True,
    widths=0.55,
    showfliers=False
)

for patch in bp["boxes"]:
    patch.set_alpha(0.6)

for median in bp["medians"]:
    median.set_linewidth(1.5)

np.random.seed(123)
for i, comp in enumerate(comparisons, start=1):
    y = all_reps[comp]
    x = np.random.normal(loc=i, scale=0.05, size=len(y))
    ax.scatter(x, y, s=28, alpha=0.8)

ax.set_ylabel("Pearson correlation")
ax.set_xlabel("")
ax.set_title("Distribution of compartment correlations across chromosomes")
ax.set_ylim(0.2, 1.02)

# 去掉横向虚线/网格线
ax.grid(False)

plt.tight_layout()
# plt.savefig("H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\matrix\\compartment\\chromosome_compartment_correlation_boxplot.png", dpi=300)
plt.savefig("/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/HiC_Rep/HUVEC_WT_LSS_OSS_reproducibility_boxplot_new_style.pdf")








    