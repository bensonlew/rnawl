# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os
import re
import argparse
import random
import pandas as pd
from pandas import Series


class calculate_pangenome(object):
    def __init__(self):
        self.all_sample_dict = {}

    def calculate(self, cluster, group, out):
        """
        根据group表生成对应的group_list
        :return:
        """
        self.sample_list = []
        with open(group, 'r') as f:#根据group表获取组内的对应关系
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    line = line.strip().split('\t')
                    sample_name = line[0]
                    self.sample_list.append(sample_name)
        data = pd.read_table(cluster, sep='\t', header=0)
        data = data.set_index("Cluster_ID")
        all_sample_index = {}
        for sample in self.sample_list:
            single_sample = data[sample]
            single_sample = single_sample[single_sample != '-']
            single_sample_index = single_sample.index
            all_sample_index[sample] = list(single_sample_index)

        with open(out, 'w') as w:
            w.write('Sample_name\tcluster_name\n')
            for sample in all_sample_index.keys():
                single_sample = all_sample_index[sample]
                w.write('{}\t{}\n'.format(sample, ';'.join(single_sample)))

if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[cluster_file]',required=True, help='input cluster file')
    parse.add_argument('-r', metavar='[group_file]', required=True, help='input group file')
    parse.add_argument('-out', metavar='[output_file]', required=True, help='input output file')
    args = parse.parse_args()
    infile = args.i
    group = args.r
    out =args.out
    m = calculate_pangenome()
    m.calculate(infile, group, out)