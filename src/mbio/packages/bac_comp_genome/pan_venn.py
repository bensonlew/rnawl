# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

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
        group_dict = {}
        group_list = []
        self.sample_list = []
        with open(group, 'r') as f:#根据group表获取组内的对应关系
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    line = line.strip().split('\t')
                    group_name = line[1]
                    sample_name = line[0]
                    if group_name not in group_list:
                        group_list.append(group_name)
                    group_dict[sample_name] = group_name
                    self.sample_list.append(sample_name) #获取group表中所有的样本名称
        new_group_dict = {}
        for gp in group_list:#获得每个group对应的样本list，方便计算单个分组内所有的unique
            group_li = []
            for key, value in group_dict.iteritems():
                if value == gp:
                    group_li.append(key)
            new_group_dict[gp] = group_li

        data = pd.read_table(cluster, sep='\t', header=0)
        data = data.set_index("Cluster_ID")
        print(group_list)
        print(new_group_dict)
        if len(group_list) == 1: ###分组个数为1，将group表中所有样本结果得到
            all_sample_index = {}
            for sample in self.sample_list:
                single_sample = data[sample]
                single_sample = single_sample[single_sample != '-']
                single_sample_index = single_sample.index
                all_sample_index[sample] = list(set(single_sample_index))

            with open(out, 'w') as w:
                w.write('Sample_name\tGene_num\tGene_name\n')
                for sample in self.sample_list:
                    single_sample = all_sample_index[sample]
                    unique_name_num = len(single_sample)
                    w.write('{}\t{}\t{}\n'.format(sample, unique_name_num, ';'.join(single_sample)))
        else:
            record_list = []
            cluster_list = []
            all_group_index = {}
            for gup in group_list:
                sample_list = new_group_dict[gup]
                for sample in sample_list:
                    single_sample = data[sample]
                    single_sample = single_sample[single_sample != '-']
                    single_sample_index = single_sample.index
                    if sample not in record_list:
                        record_list.append(sample)
                        single_sample_list = list(single_sample_index)
                        cluster_list = cluster_list + single_sample_list  #取和集
                all_group_index[gup] = list(set(cluster_list))
            with open(out, 'w') as w:
                w.write('Group_name\tGene_num\tGene_name\n')
                for grp in group_list:
                    gro_list = all_group_index[grp]
                    gro_name = ';'.join(gro_list)
                    unique_name_num = len(gro_list)
                    w.write('{}\t{}\t{}\n'.format(grp, unique_name_num, gro_name))


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



