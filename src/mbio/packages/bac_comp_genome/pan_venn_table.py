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
        根据二维关系表整理成get_homologus的输入表格
        统计出core、pan和new基因的公式所需要的数据
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
        all_sample_num = len(self.sample_list)
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
        if len(group_list) == 1:
            record_list = []
            cluster_list = []
            all_sample_index = {}

            for sample in self.sample_list:
                single_sample = data[sample]
                single_sample = single_sample[single_sample != '-']
                single_sample_index = single_sample.index
                # self.all_sample_dict[sample] = set(single_sample)
                # dict_index = {}
                # for i in single_sample_index:
                #     dict_index[i] = data[sample][i]
                # all_sample_index[sample] = dict_index
                all_sample_index[sample] = ';'.join(single_sample_index)
                if sample not in record_list:
                    record_list.append(sample)
                    single_sample_list = list(single_sample_index)
                    cluster_list = cluster_list + single_sample_list  #取和集

            dict = {}#统计出数据，得到unique总的集合
            for key in cluster_list:
                dict[key] = dict.get(key, 0) + 1
            core_list = []
            total_unique = []
            for key,value in dict.items():
                if value == 1:
                    total_unique.append(key)
                elif value == all_sample_num:
                    core_list.append(key)
            print(core_list[0:2])
            with open(out, 'w') as w:
                w.write('Sample_name\tGene_num\tGene_name\n')
                for sample in self.sample_list:
                    sample_na = sample + ' only'
                    single_sample = data[sample]
                    single_sample = single_sample[single_sample != '-']
                    single_sample_index = single_sample.index
                    unique_index = set(single_sample_index).intersection(set(total_unique))
                    # unique_name = []
                    # for index in unique_index:
                    #     index_name = all_sample_index[sample][index]
                    #     uniq_name = index_name.split("|")[1]
                    #     unique_name.append(uniq_name)
                    unique_name_num = len(unique_index)
                    w.write('{}\t{}\t{}\n'.format(sample_na, unique_name_num, ';'.join(unique_index)))
                core_list_name = " & ".join(self.sample_list)
                gene_list_name = ";".join(core_list)
                core_gene_num = len(core_list)
                w.write('{}\t{}\t{}\n'.format(core_list_name, core_gene_num, gene_list_name,))
        else:
            record_list = []
            cluster_list = []
            all_group_index = {}
            for gup in group_list:
                all_sample_index = {}
                sample_list = new_group_dict[gup]
                for sample in sample_list:
                    single_sample = data[sample]
                    single_sample = single_sample[single_sample != '-']
                    single_sample_index = single_sample.index
                    # self.all_sample_dict[sample] = set(single_sample)
                    # dict_index = {}
                    # for i in single_sample_index:
                    #     if data[sample][i] != '-':
                    #         dict_index[i] = data[sample][i]
                    # print(dict_index)
                    if sample not in record_list:
                        record_list.append(sample)
                        single_sample_list = list(single_sample_index)
                        cluster_list = cluster_list + single_sample_list  #取和集
                    all_sample_index[sample] = list(single_sample_index)
                all_group_index[gup] = all_sample_index
            # print(all_group_index["A"]['GCF_000007685.1'])
            dict = {}
            for key in cluster_list:
                dict[key] = dict.get(key, 0) + 1
            core_list = []
            total_unique = []
            for key, value in dict.items():
                if value == 1:
                    total_unique.append(key)
                elif value == all_sample_num:
                    core_list.append(key)
            with open(out, 'w') as w:
                w.write('Group_name\tGene_num\tGene_name\n')
                for grp in group_list:
                    group_na = grp + ' only'
                    sample_list = new_group_dict[grp]
                    group_index = []
                    for sample in sample_list:
                        single_sample = data[sample]
                        single_sample = single_sample[single_sample != '-']
                        single_sample_index = list(single_sample.index)
                        group_index += single_sample_index
                    # print(group_index)
                    group_uniqe_index = set(group_index)
                    unique_index = set(group_uniqe_index).intersection(set(total_unique))
                    # print(unique_index)
                    # unique_name = []
                    # for index in unique_index:
                    #     same_index = []
                    #     for sample in sample_list:
                    #         if index in all_group_index[grp][sample].keys():
                    #             index_name = all_group_index[grp][sample][index]
                    #             uniq_name = index_name.split("|")[1]
                    #             same_index.append(uniq_name)
                    #     same_index_name = ','.join(same_index)
                    #     unique_name.append(same_index_name)
                    unique_name_num = len(unique_index)
                    w.write('{}\t{}\t{}\n'.format(group_na, unique_name_num, ';'.join(unique_index)))
                core_list_name = " & ".join(group_list)
                gene_list_name = ";".join(core_list)
                core_gene_num = len(core_list)
                w.write('{}\t{}\t{}\n'.format(core_list_name, core_gene_num, gene_list_name))


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



