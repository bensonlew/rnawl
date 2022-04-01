# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import pandas as pd
import numpy as np

df = pd.DataFrame(np.random.rand(100, 10))
func = lambda row: (row - np.min(row)) / (np.max(row) - np.min(row))
zdf = pd.DataFrame(func(row) for i, row in df.iterrows())

col = zdf.columns
idx = zdf.index

mt2rs_dict = dict((c, list()) for c in col)
for i, row in zdf.iterrows():
    for k, v in row.items():
        if v == row.max():
            mt2rs_dict[k].append(row)

df_list = list()
for t in col:
    tmp_df = pd.DataFrame(mt2rs_dict[t])
    df_list.append(tmp_df.sort_values(by=col[-1]))

sdf = pd.concat(df_list)
