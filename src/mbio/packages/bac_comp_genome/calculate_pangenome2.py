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

    def calculate(self, cluster, out):
        """
        根据二维关系表整理成get_homologus的输入表格
        统计出core、pan和new基因的公式所需要的数据
        :return:
        """
        with open(cluster, 'r') as f: #获取所有样本的名称
            line = f.readline()
            line = line.strip().split("\t")
            self.sample_list = line[3:]
        pan_cluster = os.path.join(out, 'pan_cluster.xls')
        core_cluster = os.path.join(out, 'core_cluster.xls')
        new_cluster = os.path.join(out, 'new_gene_cluster.xls')
        data = pd.read_table(cluster, sep='\t', header=0)
        data1 = data.set_index("Cluster_ID")
        print(self.sample_list)
        for sample in self.sample_list:
            single_sample = data1[data1[sample] != "-"].index
            self.all_sample_dict[sample] = set(single_sample)
        outpan = open(pan_cluster, 'w')
        outcore = open(core_cluster, 'w')
        outnew = open(new_cluster, 'w')
        all_sample_num = len(self.sample_list)
        all_sample_name = []
        for i in range(1, all_sample_num+1):
            sample_name = "g" + str(i)
            all_sample_name.append(sample_name)
        outpan.write("{}\n".format("\t".join(all_sample_name)))
        outcore.write("{}\n".format("\t".join(all_sample_name)))
        outnew.write("{}\n".format("\t".join(all_sample_name)))
        for i in range(10):#获取10组数据
            core_list, pan_list, newgen_list = self.get_all_sample_data(i, all_sample_num)
            outpan.write("{}\n".format("\t".join([str(j) for j in pan_list])))
            outcore.write("{}\n".format("\t".join([str(j) for j in core_list])))
            outnew.write("{}\n".format("\t".join([str(j) for j in newgen_list])))
        outpan.close()
        outcore.close()
        outnew.close()
        # print(self.all_sample_dict)

    def get_all_sample_data(self,time, total_sample_num):
        """
        计算所有样本的需要计算的点,需要传入样本数和所有样本对应--所有cluster情况的一个字典
        :return:str(random.randint(1, 10000))
        """
        core_curve = []
        pan_curve = []
        newgene_curve = []

        for i in range(1,total_sample_num+1):
            choose_sample_list = self.get_random_new(time,i, total_sample_num)#获得需要进行计算的样本名称的list
            # core_num = self.get_cluster_list(choose_sample_list, method='intersection')
            # pan_num = self.get_cluster_list(choose_sample_list, method='union')
            newgen_num = self.get_cluster_list(choose_sample_list,method='difference')
            # core_curve.append(core_num)
            # pan_curve.append(pan_num)
            newgene_curve.append(newgen_num)
        return (core_curve,pan_curve,newgene_curve)

    def get_cluster_list(self, origin_list, method=None):
        """
        根据传入的样本list进行计数，并返回计算结果
        :return:
        """
        record_list = []
        sample_first = origin_list[0] #选取第一个作为初始的list集合
        cluster_list = self.all_sample_dict[sample_first]
        record_list.append(sample_first)
        if method in ['union', 'intersection']:
            for sample in origin_list: #递归进行选取样本
                sample_cluster_list = self.all_sample_dict[sample]
                #print(sample_cluster_list[0])
                if sample not in record_list:
                    record_list.append(sample)
                    if method in ['intersection']:
                        cluster_list = cluster_list.intersection(sample_cluster_list) #取交集
                    elif method in ['union']:
                        cluster_list = cluster_list.union(sample_cluster_list) #取并集
            cluster_num = len(cluster_list)
        else:
            cluster_list = list(cluster_list) #将集合转为list格式，然后对总的list进行操作
            for sample in origin_list: #递归进行选取样本 计算unique
                sample_cluster_list = list(self.all_sample_dict[sample])
                if sample not in record_list:
                    record_list.append(sample)
                    cluster_list = cluster_list + sample_cluster_list #取和集
            #set_cluster_list = set(cluster_list)
            cluster_num = 0
            cluster_list = Series(cluster_list)#转成series，然后对转成后的结果进行计数，随后做统计
            result = cluster_list.value_counts()
            for i in result.values:
                if i == 1:
                    cluster_num += 1
        return cluster_num

    def get_random(self, num, total_sample_num):
        """
        从随机数中获取样本名称的list
        :return:
        """
        h = set()
        while(len(h)<num): #产生随机数，得到随机的组合，从而能保证随机每次取出来的结果都不一样
            h.add(random.randint(0,total_sample_num-1))
        sample_index_list = list(h) #获取了一组随机数（这个个数与选择的样本数是相同的）
        sample_list = []
        for j in sample_index_list:
            sample_list.append(self.sample_list[j])
        return sample_list

    def get_random_new(self, time, num, total_sample_num):

        sample_index_list = [] #获取了一组随机数（这个个数与选择的样本数是相同的）
        while len(sample_index_list) < num:
            if time >= total_sample_num:
                new_time = time - total_sample_num
                sample = self.sample_list[new_time]
                sample_index_list.append(sample)
                time += 1
            else:
                sample = self.sample_list[time]
                sample_index_list.append(sample)
                time += 1
        return sample_index_list

if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[cluster_file]',required=True, help='input cluster file')
    parse.add_argument('-out', metavar='[output_dir]', required=True, help='input output dir')
    args = parse.parse_args()
    infile = args.i
    out =args.out
    m = calculate_pangenome()
    m.calculate(infile, out)



