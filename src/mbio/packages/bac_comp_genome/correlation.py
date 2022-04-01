# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import kendalltau
import numpy as np
import os
import sys


def correlation(otu_table, method):
    if not os.path.isfile(otu_table):
        raise Exception('OTU或Species文件不存在:{}'.format(otu_table))

    pvaluefile = open("{}_pvalue.xls".format(method), 'w')
    with open("correlation_{}.xls".format(method), "w") as newfile:
        df = pd.read_csv(otu_table,sep="\t", dtype=str)
        analyzed_ids = list(df.columns[1:])
        allsamples = {}
        k=0
        for samplename in analyzed_ids:
            allsamples[samplename]=k
            k+=1
        newfile.write('name'+'\t')
        pvaluefile.write('name'+'\t')
        all_factors=list(df.columns)[1:]
        for j in range(0,len(all_factors)):
            pvaluefile.write(all_factors[j])
            newfile.write(all_factors[j])
            if j!=len(all_factors)-1:
                newfile.write('\t')
                pvaluefile.write('\t')
        newfile.write('\n')
        pvaluefile.write('\n')
        for i in range(0,len(all_factors)):
            first_sample = all_factors[i]
            corr_write_line = [first_sample]
            pavalue_write_line = [first_sample]
            isnan = False
            new_values = [ float(x) for x in df[first_sample].values]
            first_array = np.asarray(new_values)
            print(first_array)
            for f in range(0,len(all_factors)):
                second_sample = all_factors[f]
                values2 = [float(x) for x in df[second_sample].values]
                second_array = np.asarray(values2)
                print(second_array)
                # correlation_result = []
                if method == "pearsonr":
                    correlation_result = pearsonr(first_array,second_array)
                if method == "spearmanr":
                    correlation_result = spearmanr(first_array,second_array)
                if method == "kendalltau":
                    correlation_result = kendalltau(first_array,second_array)
                if np.isnan(correlation_result[0]):
                    isnan = True
                corr_write_line.append(str(correlation_result[0]))
                pavalue_write_line.append(str(correlation_result[1]))
            if isnan:
                continue
            else:
                newfile.write("\t".join(corr_write_line)+"\n")
                pvaluefile.write("\t".join(pavalue_write_line)+"\n")


if __name__ == '__main__':
    correlation(sys.argv[1],sys.argv[2])
