# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from optparse import OptionParser
import pandas as pd
import os
import argparse


class Statistics(object):
    """
    此脚本用于统计不同软件的结果
    """
    def __init__(self):
        """
        初始化
        :return:
        """

    def calculate(self, input_file, method, output, step=None):
        """
        根据传入的文件夹统计结果表
        :param input_file: 文件夹路径
        :param method: 方法
        :param step: 步长
        :param output: 输出结果文件夹
        :return:
        """
        print("开始进行计算啦")
        if not os.path.exists(input_file):
            raise Exception("传入的文件夹路径不存在")
        data = pd.read_table(input_file, sep="\t", header=0)
        columns = list(data.columns)
        if "ASV ID" in columns:
            data = data.set_index("ASV ID")
        else:
            data = data.set_index("#OTU ID")
        sample_list = list(data.columns)
        n = 0 # 对样本进行计数
        with open(output, 'w') as w:
            for sample in sample_list:
                data_level = data[sample]
                data_level.columns = [sample]
                if method in ["relative"]:
                    total = data_level.sum(axis=0)
                    data_level = data_level.apply(lambda x: float(x)/total)
                data_level = data_level.sort_values(ascending=False)
                if step not in ['']:##按步长挑选数据
                    step = int(step)
                    all_list = data_level.values
                    new_list = list(all_list)

                else:
                    new_list = data_level.values
                    # step_list = range(len(new_list))
                if n == 0:## 为了写header，因为不知道抽了多少序列数，步长默认选择100
                    # index_list = [list(data_level.index)[x] for x in step_list]
                    # value_list = [list(data_level.values)[x] for x in step_list]
                    rank = range(len(new_list)) # 对选取的步长进行排序
                    w.write("sample" +"\t" +"\t".join([str(x+1) for x in rank]) + "\n") #
                    w.write(sample +"\t" +"\t".join([str(x) for x in new_list]) + "\n")
                else:
                    w.write(sample +"\t" +"\t".join([str(x) for x in new_list]) + "\n")
                n += 1
        # all_rank = pd.read_table("rank.xls", sep="\t", header=0)
        # all_rank = all_rank.T
        # all_rank.to_csv(output, index=0, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='input', type=str, help='input dir file')
    parser.add_argument('-m', metavar='method', type=str, help='calculate method')
    parser.add_argument('-s', metavar='step', type=str, help='step size')
    parser.add_argument('-o', metavar='output', type=str, help='output directory containing expression matrix table')
    args = parser.parse_args()
    input_file = args.i
    method = args.m
    output = args.o
    step = args.s
    stat = Statistics()
    try:
        stat.calculate(input_file, method,output, step=step)
    except:
        raise Exception("计算失败！")

