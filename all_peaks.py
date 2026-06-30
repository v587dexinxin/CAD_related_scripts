# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 15:55:21 2022

@author: 86182
"""

import numpy as np
import pandas as pd
import os


'''
Load files
'''

path = 'D:\\work\\lvjian\\ChIP-seq\\ChIP_seq1\\peaks\\selected_peaks\\'
all_files = [f for f in os.listdir(path)]

# all_files = ['1MLF1-273-1_peaks.broadPeak',
#  '1MLF1-27ac-3_peaks.narrowPeak',
#  '1MLF1-EZH-5_peaks.narrowPeak',
#  '1MLF1-SU-7_peaks.narrowPeak',
#  '1MLF2-273-2_peaks.broadPeak',
#  '1MLF2-27ac-4_peaks.narrowPeak',
#  '1MLF2-EZH-6_peaks.narrowPeak',
#  '1MLF2-SU-8_peaks.narrowPeak',
#  '2NC1-273-1_peaks.broadPeak',
#  '2NC1-27ac-3_peaks.narrowPeak',
#  '2NC1-EZH-5_peaks.narrowPeak',
#  '2NC1-SU-7_peaks.narrowPeak',
#  '2NC2-273-2_peaks.broadPeak',
#  '2NC2-27ac-4_peaks.narrowPeak',
#  '2NC2-EZH-6_peaks.narrowPeak',
#  '2NC2-SU-8_peaks.narrowPeak']

data_count = {}

union = pd.DataFrame()

n = 0
for k in all_files:
    file = pd.read_table(path + k , sep='\t' ,header = None)
    # file = file[~file.index.duplicated(keep='first')]
    file = file.loc[:,:2]
    file.columns = ['chr' , 'start' , 'end']
    data_count[k] = file
    union = pd.concat([union , file] , axis = 0)
    n += len(file)
    
    
union = union.sort_values(by = ['chr' , 'start'])



    
union.to_csv(path + 'all_peaks.txt' , sep = '\t' , index = False)


