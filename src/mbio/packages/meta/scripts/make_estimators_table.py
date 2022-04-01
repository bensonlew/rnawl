#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import glob



file_path = glob.glob(r"*.summary")
index_type = []
sample_name = []
line_list = []
for fp in file_path:
    sp_name = fp.split('.')[1]
    sample_name.append(sp_name)
    # print sample_name
    with open(fp, 'r') as f:
        index_type = f.readline().split('\t')
        index_type.pop(0)
        for line in f:
            line = line.strip().split('\t')
            line.pop(0)
            line = sp_name + '\t' + '\t'.join(line)
            line_list.append(line)

########################added if to else by yiru 20170421 
if os.path.exists(r"adiv_chao1_pd.txt"):
    with open("adiv_chao1_pd.txt",'r') as f:
        pdinfo = [line.rstrip('\n') for line in f]
    pddict = {}
    for i in range(1,len(pdinfo)):
        line = pdinfo[i]
        pddict[line.split('\t')[0]] = line.split('\t')[1]
    new_line_list = []
    for item in line_list:
        sp_name = item.split('\t')[0]
        item = item + '\t' + pddict[sp_name]
        new_line_list.append(item)
    index_type[-1] = index_type[-1][:-1]
    index_type.append("pd\n")
    line_list = new_line_list
else:
    pass
#########################
with open('estimators.xls', 'wb') as e:
    first_line = 'sample' + '\t' + '\t'.join(index_type)
    e.write(first_line)
    for line in line_list:
        e.write(line + '\n')
