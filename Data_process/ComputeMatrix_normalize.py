# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:45:40 2026

@author: lenovo
"""


import pandas as pd


####-------------------------HUVEC_WT_LS_OS_peaks---------------------------------

data = pd.read_csv('/scratch/2026-05-11/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/peaks_intensity_heatmap/HUVEC_union_peaks' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)

bins_num = 400

n =  len(data.columns) // bins_num

data_genes = data.iloc[:,:6]



data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    
    
    
data_new.to_csv('/scratch/2026-05-11/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/peaks_intensity_heatmap/HUVEC_union_peaks_norm' , header=None , index = None , sep = '\t')




####-------------------------HUVEC_ATAC_WT_LS_OS_peaks---------------------------------

data = pd.read_csv('/scratch/2026-05-11/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/peaks_intensity_heatmap/HUVEC_union_peaks_ATAC' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)

bins_num = 400

n =  len(data.columns) // bins_num

data_genes = data.iloc[:,:6]



data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    
    
    
data_new.to_csv('/scratch/2026-05-11/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/peaks_intensity_heatmap/HUVEC_union_peaks_ATAC_norm' , header=None , index = None , sep = '\t')



####-------------------------HUVEC_ATAC_WT_LS_OS_peaks_rescaled---------------------------------

data = pd.read_csv('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/peaks_intensity_heatmap/HUVEC_union_peaks_ATAC_rescaled' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)

bins_num = 400

n =  len(data.columns) // bins_num

data_genes = data.iloc[:,:6]



data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    
    
    
data_new.to_csv('/scratch/2026-05-18/bio-shenw/Cardiovascular_disease_STARR-seq/HUVEC_Cardiovascular_disease_moudle/HiRPC/plots/peaks_intensity_heatmap/HUVEC_union_peaks_ATAC_rescaled_norm' , header=None , index = None , sep = '\t')







