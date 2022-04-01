# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from __future__ import print_function

import logging
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('usage: {} <bam_depth> <range> <sample_name> <output_dir>'.format(os.path.basename(sys.argv[0])))
    bam_depth = sys.argv[1]
    range = int(sys.argv[2])
    sample = sys.argv[3]
    output_dir = sys.argv[4]
    df = pd.read_table(bam_depth, header=None)
    max = df[1].max()
    a = range
    data = dict()
    while a <= df[1].max():
        data[a] = list(df[(df[1] <= a) & (df[1] > a-range)][2])
        a += range
    # if df[1].max() > a-range:
    #     col = ">" + str(a-range)
    #     data[col] = list(df[df[1] >= a][2])
    plot_data = pd.DataFrame.from_dict(data)
    plot_data_sort = pd.DataFrame(plot_data, columns=sorted(data.keys()))
    plt.figure(figsize=(9, 9))
    plot_data_sort.boxplot(sym='ro', patch_artist=True, rot=45)
    title = sample + " coverage"
    plt.title(title, fontsize=20, fontweight='bold')
    plt.xlabel("bp bins", fontsize=14, fontweight='bold')
    plt.ylabel("# of reads aligned", fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(output_dir, sample + '_coverage.png'))
    plt.savefig(os.path.join(output_dir, sample + '_coverage.pdf'))
