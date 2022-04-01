# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd


class calculate_pangenome(object):
    def __init__(self):
        self.all_sample_dict = {}

    def func(self, x, a, b):
        return a*x**b

    def calculate(self, infile, outfile, method):
        """
        根据二维关系表整理成get_homologus的输入表格
        统计出core、pan和new基因的公式所需要的数据
        :return:
        """
        data = pd.read_table(infile, sep='\t', header=0)
        columns = data.columns
        xdata = [x for x in range(1, len(columns)+1)]
        ydata = []
        for sample in columns:
            if method in ["median"]:
                sample_method = np.median(data[sample])
            elif method in ['average']:
                sample_method = np.mean(data[sample])
            elif method in ["bincount"]:#众数
                counts = np.bincount(data[sample])
                sample_method = np.argmax(counts)
            else:
                sample_method = data[sample][0] #取样本选择的第一个数据
            ydata.append(sample_method)
        # ydata = ydata[1:]
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        print(xdata)
        print(ydata)
        popt, pcov = curve_fit(self.func, xdata, ydata)
        print(popt)
        print(pcov)
        outf = open(outfile, 'w')
        outf.write("# newgene fit converged\n")
        outf.write("# Formula: \n ~ newgenes(g) = {} * g ** {} \n".format(str(popt[0]), str(popt[1])))
        outf.write("# fitted values (genomes, genes) :\n")
        new_ydata = []
        # for x in xdata:
        #     y = popt[0] * x ** popt[1]
        #     new_ydata.append(y)
        aa = popt[0]
        bb = popt[1]
        new_ydata = self.func(xdata, aa, bb)
        print(ydata)
        print(new_ydata)
        outf.write("{}\n".format('\t'.join(str(x) for x in xdata)))
        # outf.write("{}\n".format('\t'.join(str("%.3f" %(y)) for y in new_ydata)))
        outf.write("{}\n".format('\t'.join(str(y) for y in new_ydata)))
        outf.close()


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[cluster_file]',required=True, help='input cluster file')
    parse.add_argument('-o', metavar='[output_file]', required=True, help='input output_file')
    parse.add_argument('-m', metavar='[method]', required=True, help='input method')
    args = parse.parse_args()
    infile = args.i
    out = args.o
    method = args.m
    m = calculate_pangenome()
    m.calculate(infile, out, method)



