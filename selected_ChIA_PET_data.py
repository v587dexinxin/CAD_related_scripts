# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 16:17:10 2023

@author: 86182
"""

import numpy as np
from itertools import islice
import pandas as pd



def write_tofiles(data , outfile):
    '''
    '''
    out = open(outfile , 'w')
    # out.writelines('\t'.join([str(x) for x in data.columns]))
    for i in data.index:
        out.writelines('\t'.join([str(x) for x in data.loc[i]]) + '\n')
    out.close()
    

data = pd.read_table('D:/work/literature_data/HCT116/ChIA-PET_Encode/LHH0127V.e500.clusters.cis__HCT116_chiapet_rep2.BE3' , header = None)
data.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF']


E2_1K = (22078285,22079431)
E3_1K= (22084459,22085712)
E5_2K= (22102429,22104539)       
E6_1K= (22106931,22108246)
E7_1K= (22111222,22112368)
E8_2K= (22118010,22119777)
E9_1K= (22127135,22128465)

intervals = [E2_1K, E3_1K , E5_2K , E6_1K , E7_1K , E8_2K , E9_1K]
d_type = np.dtype({'names':['start' , 'end'] , 'formats':[np.int64 , np.int64]})
intervals = np.array(intervals, dtype = d_type)


tmp_data = data[(data['chr1'] == 'chr9') & (data['chr2'] == 'chr9')]


# selected_data = []
# for i in tmp_data.index:
#     chr1 = tmp_data.loc[i].chr1
#     chr2 = tmp_data.loc[i].chr2
#     start1 = tmp_data.loc[i].start1
#     start2 = tmp_data.loc[i].start2
#     end1 = tmp_data.loc[i].end1
#     end2 = tmp_data.loc[i].end2
#     IF = tmp_data.loc[i].IF
    
#     mask = ((end1 >= intervals['start']) & (start1 <= intervals['end']) ) | ((end2 >= intervals['start']) & (start2 <= intervals['end']) )
#     if mask.sum() > 0:
#         selected_data.append((chr1 , start1 , end1 , chr2 , start2 , end2 , IF))
            




            
selected_data = pd.DataFrame(columns=tmp_data.columns)         
for i in intervals:
    start = i[0]
    end = i[1]
    mask = ((tmp_data['start1'] <= end) & (tmp_data['end1'] >= start)) | ((tmp_data['start2'] <= end) & (tmp_data['end2'] >= start))
    if mask.sum() > 0:
        tmp = tmp_data[mask]
        selected_data = selected_data.append(tmp)
        
        
selected_data.drop_duplicates(keep = 'first' , inplace=True)
    
write_tofiles(selected_data , 'D:\work\literature_data\HCT116\ChIA-PET_Encode\\Selected_9p21_HCT116_chiapet_rep2.bedpe')

    
        
















