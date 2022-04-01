# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180925


import sys
import argparse
import os
import pandas as pd
import shutil
import scipy
from scipy import stats


class Cal_Corretation_fd(object):
    def __init__(self):
        self.sel_length_table = ""

    def correlation(self, profile1, profile2, coefficient, c_cut, p_cut, outfile, trans):
        profile_table1 = pd.read_table(profile1, sep='\t', header=0, index_col=0)
        if trans:
            profile_table2 = pd.read_table(profile2, sep='\t', header=0, index_col=0)
            print profile_table2.head()
            profile_table2 = profile_table2.T
            profile_table2.columns = profile_table2.columns.astype("string")
        else:
            profile_table2 = pd.read_table(profile2, sep='\t', header=0, index_col=0)
        # table_matrix = table.T.corr(coefficient)
        length1 = len(profile_table1)
        length2 = len(profile_table2)
        node1 = []
        node2 = []
        correlations = []
        pvalues = []
        common_sam = profile_table1.columns & profile_table2.columns
        profile_table1 = profile_table1[common_sam]
        profile_table2 = profile_table2[common_sam]
        '''
        if length < 2:
            raise Exception('丰度表中物种或功能小于2，无法计算相关性!')
        '''
        for i in range(0, length1 - 1):
            for j in range(0, length2):
                (correlation, pvalue) = self.coeff_cal(coefficient, profile_table1.iloc[i], profile_table2.iloc[j])
                node1.append(profile_table1.index[i])
                node2.append(profile_table2.index[j])
                correlations.append(correlation)
                pvalues.append(pvalue)
        tmp_frame = {"Node1": node1, "Node2": node2, "Correlation": correlations, "P_value": pvalues}
        frame = pd.DataFrame(tmp_frame, columns=["Node1", "Node2", "Correlation", "P_value"])
        frame = frame[(abs(frame["Correlation"]) >= float(c_cut)) & (frame["P_value"] <= float(p_cut))]
        frame.to_csv(outfile, sep="\t", index=False)

    def coeff_cal(self, coefficient, x, y):
        result = []
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
    parser.add_argument('-i1', metavar='[profile file1]', required=True, help='Input profile file')
    parser.add_argument('-i2', metavar='[profile file2]', required=True, help='Input profile file2')
    parser.add_argument('-o', metavar='[outputfile]', required=True, help='output Directory name')
    parser.add_argument('-c', metavar='coefficient', help='spearman, pearson or kendall', default="spearman")
    parser.add_argument('-p_cut', metavar='p-value', help='input p-value cutoff', default=0.05)
    parser.add_argument('-c_cut', metavar='coefficient value', help='input coefficient value cutoff', default=0.6)
    parser.add_argument('--trans', action='store_true', help="transform table2 to calculate")
    args = parser.parse_args()
    trans = False
    profile1 = args.i1
    profile2 = args.i2
    outfile = args.o
    coefficient = args.c
    c_cut = args.c_cut
    p_cut = args.p_cut
    if args.trans:
        trans = True
    run_cal = Cal_Corretation_fd()
    run_cal.correlation(profile1, profile2, coefficient, c_cut, p_cut, outfile, trans)
