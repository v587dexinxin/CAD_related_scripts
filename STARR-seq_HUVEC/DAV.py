# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:44:39 2026

@author: lenovo
"""

import pandas as pd
import numpy as np



HUVEC_nonrisk_e_speci = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\nonrisk_VS_risk_specific_RE\\HUVEC_nonrisk_enhancer_speci.csv' , header = 0)
HUVEC_nonrisk_s_speci = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\nonrisk_VS_risk_specific_RE\\HUVEC_nonrisk_silencer_speci.csv' , header = 0)
HUVEC_risk_e_speci = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\nonrisk_VS_risk_specific_RE\\HUVEC_risk_enhancer_speci.csv' , header = 0)
HUVEC_risk_s_speci = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\STARR-seq\\DiffBind\\Regulatory_Element\\edgR\\nonrisk_VS_risk_specific_RE\\HUVEC_risk_silencer_speci.csv' , header = 0)



nonrisk_e_speci_seq = set(HUVEC_nonrisk_e_speci['seq'])
nonrisk_s_speci_seq = set(HUVEC_nonrisk_s_speci['seq'])
risk_e_speci_seq = set(HUVEC_risk_e_speci['seq'])
risk_s_speci_seq = set(HUVEC_risk_s_speci['seq'])


len(nonrisk_e_speci_seq | nonrisk_s_speci_seq | risk_e_speci_seq | risk_s_speci_seq)

