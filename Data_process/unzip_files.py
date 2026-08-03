# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 17:09:48 2022

@author: 86182
"""


import gzip
import os



def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "")
    #获取文件的名称，去掉
    g_file = gzip.GzipFile(file_name)
    a = g_file.read()
    s = str(a, encoding="utf-8")
    s = s.replace('\r', '')
    out = open(f_name, "w")
    out.writelines(s)
    g_file.close()
    out.close()

    

file_path = 'D:\\work\\Postdoctoral\\GWAS疾病位点检测\\literature\\hapmap\\'    
file_names = os.listdir(file_path)

for f in file_names:
    un_gz(os.path.join(file_path , f))
    

    
    
    




