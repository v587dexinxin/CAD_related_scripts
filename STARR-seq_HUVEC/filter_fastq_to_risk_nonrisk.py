# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 16:32:33 2025

@author: lenovo
"""

import pandas as pd
import numpy as np
import os
import gzip
from Bio import SeqIO



#prefixes=("ATCACG", "CGATGT", "NTCACG", "NGATGT")
risk_pref = ("CGATGT", "NGATGT")
nonrisk_pref = ("ATCACG", "NTCACG")



def filter_paired_fastq(input_r1, input_r2, input_r3, r_output_r1, r_output_r2, r_output_r3, n_output_r1, n_output_r2, n_output_r3, risk_pref = ("CGATGT", "NGATGT") , nonrisk_pref = ("ATCACG", "NTCACG")):
    with gzip.open(input_r1, "rt") as infile_r1, gzip.open(input_r2, "rt") as infile_r2, gzip.open(input_r3, "rt") as infile_r3,\
         open(r_output_r1, "wt") as r_outfile_r1, open(r_output_r2, "wt") as r_outfile_r2 , open(r_output_r3, "wt") as r_outfile_r3, \
         open(n_output_r1, "wt") as n_outfile_r1, open(n_output_r2, "wt") as n_outfile_r2 , open(n_output_r3, "wt") as n_outfile_r3:
        
        records_r1 = SeqIO.parse(infile_r1, "fastq")
        records_r2 = SeqIO.parse(infile_r2, "fastq")
        records_r3 = SeqIO.parse(infile_r3, "fastq")
        
        
        
        for rec1, rec2, rec3 in zip(records_r1, records_r2, records_r3):
            if rec1.seq.startswith(risk_pref):
                SeqIO.write(rec1 , r_outfile_r1 , 'fastq')
                SeqIO.write(rec2 , r_outfile_r2 , 'fastq')
                SeqIO.write(rec3 , r_outfile_r3 , 'fastq')


            elif rec1.seq.startswith(nonrisk_pref):
                SeqIO.write(rec1 , n_outfile_r1 , 'fastq')
                SeqIO.write(rec2 , n_outfile_r2 , 'fastq')
                SeqIO.write(rec3 , n_outfile_r3 , 'fastq')
                
                
            else:
                pass
            
        r_outfile_r1.close()
        r_outfile_r2.close()
        r_outfile_r3.close()
        
        n_outfile_r1.close()
        n_outfile_r2.close()
        n_outfile_r3.close()
        
        
                
        
        
        

for sample in ['HCT116-R1-4_S1_L003' , 'HCT116-P-R2-2_S6_L003' , 'HCT116-R2-8_S2_L003']:
    input_r1 = "/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R1_001.fastq.gz"
    input_r2 = "/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R2_001.fastq.gz"
    input_r3 = "/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R3_001.fastq.gz"
    
    r_output_r1 = input_r1.split("/")[-1].split('_001')[0] + '_risk' + input_r1.split("/")[-1].split('_001')[1].rstrip('.gz')
    r_output_r2 = input_r2.split("/")[-1].split('_001')[0] + '_risk' + input_r2.split("/")[-1].split('_001')[1].rstrip('.gz')
    r_output_r3 = input_r3.split("/")[-1].split('_001')[0] + '_risk' + input_r3.split("/")[-1].split('_001')[1].rstrip('.gz')
    
    n_output_r1 = input_r1.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r1.split("/")[-1].split('_001')[1].rstrip('.gz')
    n_output_r2 = input_r2.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r2.split("/")[-1].split('_001')[1].rstrip('.gz')
    n_output_r3 = input_r3.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r3.split("/")[-1].split('_001')[1].rstrip('.gz')
    filter_paired_fastq(input_r1, input_r2, input_r3, r_output_r1, r_output_r2, r_output_r3, n_output_r1, n_output_r2, n_output_r3, risk_pref = ("CGATGT", "NGATGT") , nonrisk_pref = ("ATCACG", "NTCACG"))
    
    
    
    
for sample in ['HUVEC-P-R1-5_S7_L003' , 'HUVEC-P-R2-7_S8_L003' , 'HUVEC-R1-3_S3_L003' , 'HUVEC-R2-6_S4_L003']:
    input_r1 = "/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R1_001.fastq.gz"
    input_r2 = "/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R2_001.fastq.gz"
    input_r3 = "/scratch/2025-03-24/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R3_001.fastq.gz"
    
    r_output_r1 = input_r1.split("/")[-1].split('_001')[0] + '_risk' + input_r1.split("/")[-1].split('_001')[1].rstrip('.gz')
    r_output_r2 = input_r2.split("/")[-1].split('_001')[0] + '_risk' + input_r2.split("/")[-1].split('_001')[1].rstrip('.gz')
    r_output_r3 = input_r3.split("/")[-1].split('_001')[0] + '_risk' + input_r3.split("/")[-1].split('_001')[1].rstrip('.gz')
    
    n_output_r1 = input_r1.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r1.split("/")[-1].split('_001')[1].rstrip('.gz')
    n_output_r2 = input_r2.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r2.split("/")[-1].split('_001')[1].rstrip('.gz')
    n_output_r3 = input_r3.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r3.split("/")[-1].split('_001')[1].rstrip('.gz')
    filter_paired_fastq(input_r1, input_r2, input_r3, r_output_r1, r_output_r2, r_output_r3, n_output_r1, n_output_r2, n_output_r3, risk_pref = ("CGATGT", "NGATGT") , nonrisk_pref = ("ATCACG", "NTCACG"))
    
    






    
for sample in ['HCT116-P-R2-2_S6_L003']:
    input_r1 = "/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R1_001.fastq.gz"
    input_r2 = "/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R2_001.fastq.gz"
    input_r3 = "/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/data/231218_E00517_0917_BHJV33CCX2/fastq/" + sample + "_R3_001.fastq.gz"
    
    r_output_r1 = input_r1.split("/")[-1].split('_001')[0] + '_risk' + input_r1.split("/")[-1].split('_001')[1].rstrip('.gz')
    r_output_r2 = input_r2.split("/")[-1].split('_001')[0] + '_risk' + input_r2.split("/")[-1].split('_001')[1].rstrip('.gz')
    r_output_r3 = input_r3.split("/")[-1].split('_001')[0] + '_risk' + input_r3.split("/")[-1].split('_001')[1].rstrip('.gz')
    
    n_output_r1 = input_r1.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r1.split("/")[-1].split('_001')[1].rstrip('.gz')
    n_output_r2 = input_r2.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r2.split("/")[-1].split('_001')[1].rstrip('.gz')
    n_output_r3 = input_r3.split("/")[-1].split('_001')[0] + '_nonrisk' + input_r3.split("/")[-1].split('_001')[1].rstrip('.gz')
    filter_paired_fastq(input_r1, input_r2, input_r3, r_output_r1, r_output_r2, r_output_r3, n_output_r1, n_output_r2, n_output_r3, risk_pref = ("CGATGT", "NGATGT") , nonrisk_pref = ("ATCACG", "NTCACG"))
            
        

























