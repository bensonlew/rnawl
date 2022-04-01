# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd

"""
此脚本用于合并不同样本的丰度表
"""

input_dir = sys.argv[1]
output_path = sys.argv[2]
out_asv_table = os.path.join(output_path, "ASV_table.xls")
listdirs = os.listdir(input_dir)
df_empty = pd.DataFrame(columns=["#OTU ID"])
for file in listdirs:
    file_path = os.path.join(input_dir, file)
    line_num = 0
    with open(file_path, 'rb') as f:
        for line in f:
            if line.startswith('#OTU ID') or line.startswith('ASV ID'):
                break
            line_num += 1
    f.close()
    df = pd.read_table(file_path, sep="\t", header=line_num)
    df.rename(columns={list(df)[0]:'#OTU ID'}, inplace = True)
    data = pd.merge(df_empty, df, left_on="#OTU ID", right_on="#OTU ID", how="outer")
    df_empty = data
#df_empty.fillna("0.0", inplace=True)
df_empty.to_csv(out_asv_table, sep="\t", na_rep="0.0", index=0)