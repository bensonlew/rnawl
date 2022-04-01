# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# last_modify:20181115

import pandas as pd
import argparse
import os
from mbio.packages.statistical.metastat import two_sample_test


class SecStat(object):
    def __init__(self):
        super(SecStat, self).__init__()
        self.summary_table_name = "summary.txt"
        self.fisher_input_name = "fisher.input"
        self.fisher_output_name = "fisher.output"
        self.summary_head_list = ['total','true','false','true_percent']  # 修改此變量來改變表格列的順序
        self.gene_set = set()  # 包含的預測陽性結果的

    def run(self, anno_table, prof_table, level, out_dir, mul_test, gene_count):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        out_summary = os.path.join(out_dir, self.summary_table_name)
        fisher_input_file = os.path.join(out_dir, self.fisher_input_name)
        fisher_output_file = os.path.join(out_dir, self.fisher_output_name)
        anno = pd.read_table(anno_table, index_col=0)
        prof = pd.read_table(prof_table, index_col=0)
        if gene_count:  # 将gene的丰度转换成个数
            prof["Total"] = 0
            prof[prof > 0] = 1
            prof["Total"] = prof.sum(1)
            print(prof.head())
        summary = prof.groupby(self.filter_sec).sum().T
        if "Total" in summary.index:
            summary.drop(["Total"], inplace=True)
        summary["total"] = summary['false'] + summary['true']
        summary['true_percent'] = summary['true'] / summary["total"] * 100
        summary[self.summary_head_list].to_csv(out_summary, sep="\t", index=True, float_format='%.4f')
        anno_data = pd.concat([anno[level], prof["Total"]], axis=1).dropna()
        fisher_data = anno_data.groupby([level, self.filter_sec]).sum().unstack()["Total"].fillna(0).astype("int")
        fisher_data['total'] = fisher_data['true'] + fisher_data['false']
        tmp = fisher_data[['true', 'total']]
        csv_data = tmp[tmp["true"] != 0].copy()
        csv_data.loc["no_sec"] = tmp[tmp["true"] == 0].sum()
        csv_data.to_csv(fisher_input_file, sep="\t", index=True)
        try:
            two_sample_test(fisher_input_file,fisher_output_file,"fisher", "true", "total", mul_test=mul_test)
        except Exception,e:
            raise Exception(e)

    def add_sec_table(self, predict_table):
        pred = pd.read_table(predict_table, index_col=0).index
        self.gene_set.update(pred)

    def filter_sec(self, gene_id):
        if gene_id in self.gene_set:
            return "true"
        else:
            return "false"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', metavar='[profile file]', required=True, help='Input profile file')
    parser.add_argument('-t', metavar='[taxon file]', required=True, help='Input nr taxon file')
    parser.add_argument('-s', metavar='[sec predict file]', required=True, help='Input sec predict files, seperate by comma')
    parser.add_argument('-l', metavar='[taxon level selected]', help='choose taxon level in taxon file', default="Species")
    parser.add_argument('-m', metavar='[multi test method]', help='how to calculate corrected p value', default="fdr")
    parser.add_argument('-o', metavar='[output file]', required=True, help='output Directory name')
    parser.add_argument('-e', action='store_true', help='for gene num not gene profile')
    args = parser.parse_args()
    obj = SecStat()
    for predict_table in args.s.split(","):
        obj.add_sec_table(predict_table)
    obj.run(args.t, args.p, args.l, args.o, args.m, args.e)
