# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import os
import sys
import pandas as pd


if len(sys.argv) != 4:
    exit('USAGE: %s list infile out' %sys.argv[0])

list = sys.argv[1]
infile = sys.argv[2]
out = sys.argv[3]

with open(list, 'r') as l_r:
    # sample_list = l_r.read().split('\n')
    sample_list = [x.strip() for x in l_r.readlines()]
    while '' in sample_list:
        sample_list.remove('')

df = pd.read_table(infile, header =0, sep = '\t', index_col=0)
sort_df = df[sample_list]
sort_df.index.name = df.index.name
sort_df.to_csv(out, sep = '\t',index=True)
