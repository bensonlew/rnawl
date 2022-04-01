# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import argparse
import pandas as pd
import numpy as np
import os


class Statistics(object):
    """
    此脚本用于统计不同软件的结果
    """
    def __init__(self):
        self.mehod = {
            'tpm': {"salmon": "TPM", "kallisto": "tpm"},
            'fpkm': "FPKM"
        }
        self.file_names = {
            "kallisto": ".abundance.tsv",
            "salmon": ".quant.sf"
        }
        self.counts_name = {
            "kallisto": "est_counts",
            "salmon": "NumReads"
        }

    def calculate(self, dirs, software, method, output):
        """
        根据传入的文件夹统计结果表
        :param dirs: 文件夹路径
        :param software: 软件名称
        :param method: 方法
        :param output: 输出结果文件夹
        :return:
        """
        print("开始进行计算啦")
        if not os.path.exists(dirs):
            raise Exception("传入的文件夹路径不存在")
        all_gene_list = set()
        sample_list = []
        gene_sample_dict = {}
        gene_count_dict = {}
        files = os.listdir(dirs)
        for file in files:
            sample_path = os.path.join(dirs, file)
            my_method = self.mehod[method][software]
            file_name = self.file_names[software]
            sample_name = file.strip(file_name)
            if sample_name not in sample_list:
                sample_list.append(sample_name)
            data = pd.read_table(sample_path, sep='\t', header=0)
            method = data[my_method]
            count_name = self.counts_name[software]
            if software in ['kallisto']:
                gene_list = list(data['target_id'])
                result_dic = data.set_index('target_id').to_dict()[my_method]
                count_dic = data.set_index('target_id').to_dict()[count_name]
            else:
                gene_list = list(data['Name'])
                result_dic = data.set_index('Name').to_dict()[my_method]
                count_dic = data.set_index('Name').to_dict()[count_name]
            gene_sample_dict[sample_name] = result_dic
            gene_count_dict[sample_name] = count_dic
            all_gene_list = set(all_gene_list.union(set(gene_list)))
        out_file = os.path.join(output, "gene.tpm.matrix.xls")
        print(sample_list)
        with open(out_file, 'w') as outw:
            outw.write("Gene_id\t{}\n".format('\t'.join(sample_list)))
            for gene in all_gene_list:
                each_gene = list()
                for sample in sample_list:
                    if sample in gene_sample_dict.keys():
                        gene_dict = gene_sample_dict[sample]
                        if gene in gene_dict.keys():
                            gene_num = str(gene_dict[gene])
                            each_gene.append(gene_num)
                        else:
                            each_gene.append(str(0.000))
                gene_string = "\t".join(each_gene)
                outw.write(gene + "\t" +gene_string + "\n")
        counts_file = os.path.join(output, "gene.counts.matrix.xls")
        with open(counts_file, 'w') as outm:
            outm.write("Gene_id\t{}\n".format('\t'.join(sample_list)))
            for gene in all_gene_list:
                each_count = list()
                for sample in sample_list:
                    if sample in gene_sample_dict.keys():
                        count_dict = gene_count_dict[sample]
                        if gene in count_dict.keys():
                            count_num = str(count_dict[gene])
                            each_count.append(count_num)
                        else:
                            each_count.append(str(0.000))
                gene_counts_string = "\t".join(each_count)
                outm.write(gene + "\t" +gene_counts_string + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='input', required=True, help='input dir file')
    parser.add_argument('-m', metavar='method', required=True, help='quantitative method [tpm, fpkm]')
    parser.add_argument('-s', metavar='software', required=True, help='quantitative method [salmon, kallisto]')
    parser.add_argument('-o', metavar='output', required=True, help='output directory containing expression matrix table')
    args = parser.parse_args()

    input_dir = args.i
    method = args.m
    software = args.s
    output = args.o
    stat = Statistics()
    try:
        stat.calculate(input_dir, software, method, output)
    except:
        raise Exception("合并计算失败！")

