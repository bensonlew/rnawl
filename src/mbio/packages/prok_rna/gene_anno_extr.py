# -*- coding: utf-8 -*-
"""
@time    : 2018/10/25 9:58
@file    : gene_anno_extr.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
import argparse
import os
import sys

import pandas as pd


class GeneAnnoExtrParam(object):
    def __init__(self):
        self.outdir = None
        self.is_ready = False
        self.__parser = argparse.ArgumentParser(
            __doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    def _args_check(self, args_obj):
        self.gene_list = args_obj.gene_list
        if self.gene_list is None or not os.path.isfile(args_obj.gene_list):
            self.is_ready = False

        self.annotation_file = args_obj.annotation_file
        if self.annotation_file is None or not os.path.isfile(args_obj.annotation_file):
            self.is_ready = False

        self.outdir = args_obj.outdir
        self.filename = args_obj.outfile_name
        self.outfile = os.path.join(self.outdir, self.filename)

        self.is_ready = True

    def _add_params(self):
        self.__parser.add_argument('-g', '--gene_list', type=str, required=True,
                                   help='gene list file which has one field')
        self.__parser.add_argument('-a', '--annotation_file', type=str, required=True,
                                   help='annotation text file having header')
        self.__parser.add_argument('-f', '--outfile_name', required=True)
        self.__parser.add_argument('-o', '--outdir', help='output directory', required=True,)

    def parse_args(self, args_num):
        self._add_params()

        argvs = [e.lower() for e in sys.argv]
        if '-h' in argvs or '--help' in argvs or len(argvs) < args_num:
            self.__parser.parse_args(['-h'])

        args_obj = self.__parser.parse_args()
        self._args_check(args_obj)


class GeneAnnoExtrParamPackage(object):
    def __init__(self, param_obj):
        self.__param_obj = param_obj

    def __extr_gene(self):
        gene_set = set()
        with open(self.__param_obj.gene_list, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                gene_set.add(line)
        return gene_set

    def __get_anno(self, gene_list):
        df = pd.read_table(
            self.__param_obj.annotation_file, sep='\t', header=0, index_col=0)
        df = df.dropna(how='all').fillna('-')
        comm_genes = set(df.index) & gene_list
        df = df.loc[comm_genes, :]
        df.to_csv(self.__param_obj.outfile, sep='\t', header=True, index=True)


    def run(self):
        gene_set = self.__extr_gene()
        self.__get_anno(gene_set)


if __name__ == '__main__':
    param_obj = GeneAnnoExtrParam()
    param_obj.parse_args(4)
    if param_obj.is_ready:
        package_tool = GeneAnnoExtrParamPackage(param_obj=param_obj)
        package_tool.run()
