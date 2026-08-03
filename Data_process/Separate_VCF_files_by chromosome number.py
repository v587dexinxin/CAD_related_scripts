# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:57:46 2023

@author: 86182
"""

import pandas as pd
import numpy as np
from itertools import islice




##提取VCF中的有用信息，分别写入和按照染色体号分开写入文件中

data = open('/scratch/2023-05-08/bio-shenw/Lixinxin/ref/Human/SNP/GCF_000001405.40.vcf' , 'r')
out_all = open('/scratch/2023-05-08/bio-shenw/Lixinxin/ref/Human/SNP/GCF_000001405.40_new.vcf' , 'w') 
out_all.writelines('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT']) + '\n')

data_new = {} ; count = 0
for i in islice(data , 38 , None):
    count += 1
    if count % 1000000 == 0:
        print(count)
    i = i.strip().split('\t')
    out_all.writelines('\t'.join([i[0] , i[1] , i[2] , i[3] , i[4]]) + '\n')
    if i[0] not in data_new.keys():    
        try:
            out.close()
        except:
            pass
        data_new[i[0]] = 0
        out = open('/scratch/2023-05-08/bio-shenw/Lixinxin/ref/Human/SNP/GCF_000001405.40_Chr' + i[0] + '.vcf' , 'w')
        out.writelines('\t'.join([i[0] , i[1] , i[2] , i[3] , i[4]]) + '\n')
    else:
        out.writelines('\t'.join([i[0] , i[1] , i[2] , i[3] , i[4]]) + '\n')

out.close()
out_all.close()






    
    
    