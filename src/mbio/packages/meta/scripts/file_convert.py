# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# __version__ = 'v1.0'
# __last_modified__ = '20191122'


import pandas as pd
import sys


def remove_0_otu(otu_table):
    """
    根据传入的otu表过滤掉所有样本总和为0
    """
    data = pd.read_table(otu_table, sep='\t', header=0)
    data_samples = []
    for i in data.columns:
        if i in ["OTU ID", "KO ID", "K_desc", "CategoryL3", "PATHWAY", "Level2", "Level1", "CategoryL2",
                 "CategoryL1", "COG ID", "#Category", "Description"]:
            pass
        else:
            if i not in data_samples:
                data_samples.append(i)
    temp = data[data_samples]
    data['sum'] = temp.sum(axis=1)

    data = data[(data['sum']) != 0]
    del data['sum']
    data.to_csv(otu_table, sep='\t', index=0)


if __name__ == "__main__":
    args = sys.argv[1:]
    otu_file = args[0]
    remove_0_otu(otu_file)