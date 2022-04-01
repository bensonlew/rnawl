# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20171223


import sys
import argparse
import os
import pandas as pd
import shutil
import scipy
from scipy import stats


class Cal_Corretation(object):
    def __init__(self):
        self.sel_length_table = ""

    def correlation(self, profile, coefficient, c_cut, p_cut, outfile):
        profile = pd.read_table(profile, sep='\t', header=0, index_col=0)
        # table_matrix = table.T.corr(coefficient)
        length = len(profile)
        node1 = []
        node2 = []
        correlations = []
        pvalues = []
        if length < 2:
            raise Exception('丰度表中物种或功能小于2，无法计算相关性!')
        for i in range(0, length - 1):
            for j in range(i, length - 1):
                (correlation, pvalue) = self.coeff_cal(coefficient, profile.iloc[i], profile.iloc[j + 1])
                node1.append(profile.index[i])
                node2.append(profile.index[j + 1])
                correlations.append(correlation)
                pvalues.append(pvalue)
        tmp_frame = {"Node1": node1, "Node2": node2, "Correlation": correlations, "P_value": pvalues}
        frame = pd.DataFrame(tmp_frame, columns=["Node1", "Node2", "Correlation", "P_value"])
        frame = frame[(abs(frame["Correlation"]) >= float(c_cut)) & (frame["P_value"] <= float(p_cut))]
        frame.to_csv(outfile, sep="\t", index=False)

    def coeff_cal(self, coefficient, x, y):
        result=[]
        if coefficient == "spearman":
            result = scipy.stats.spearmanr(x, y)
        elif coefficient == "pearson":
            result = scipy.stats.pearsonr(x, y)
        elif coefficient == "kendall":
            result = scipy.stats.kendalltau(x, y)
        correlation = result[0]
        pvalue = result[1]
        return (correlation, pvalue)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('-gl', metavar='[gene_length_table]', required=True, help='Input gene_length_table')
    parser.add_argument('-i', metavar='[profile file]', required=True, help='Input profile file')
    #parser.add_argument('-rp', metavar='[gene_relative_profile]', required=True, help='Input gene_relative_profile')
    parser.add_argument('-o', metavar='[outputfile]', required=True, help='output Directory name')
    #parser.add_argument('-s', metavar='[select_genes_file]', required=True, help='Input select genes file')
    parser.add_argument('-c', metavar='coefficient', help='spearman, pearson or kendall', default="spearman")
    parser.add_argument('-p_cut', metavar='p-value', help='input p-value cutoff', default=0.05)
    parser.add_argument('-c_cut', metavar='coefficient value', help='input coefficient value cutoff', default=0.6)
    args = parser.parse_args()
    profile = args.i
    outfile = args.o
    coefficient = args.c
    c_cut = args.c_cut
    p_cut = args.p_cut
    run_cal = Cal_Corretation()
    run_cal.correlation(profile, coefficient, c_cut, p_cut, outfile)
