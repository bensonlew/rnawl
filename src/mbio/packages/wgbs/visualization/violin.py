# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from __future__ import print_function

import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: {} <result FILE> <output DIR>'.format(os.path.basename(sys.argv[0])))
    result_fp = sys.argv[1]
    output_dir = sys.argv[2]
    df = pd.read_table(result_fp)
    if 'context' in df:
        df = df[df['context'] == 'CG']
    ratio_columns = [col_name for col_name in df if col_name.endswith('_ratio')]
    ratio_columns.sort()
    ratio_df = df.reindex(ratio_columns, axis=1)
    ratio_df.rename({col_name: col_name[:-6].split('-')[2] for col_name in ratio_columns}, axis=1, inplace=True)
    ratio_df.fillna(0, inplace=True)
    if 'context' not in df:
        ratio_df /= 100
    data = list()
    for col_name in ratio_df:
        col = ratio_df[col_name]
        data.append(col[col != 0].tolist())
    _, ax = plt.subplots()
    ax.violinplot(data)
    ax.set_ylabel('methylation ratio')
    ax.set_title('Methylation distribution')
    plt.setp(ax, xticks=np.arange(ratio_df.shape[1]) + 1, xticklabels=list(ratio_df.columns))
    for xtick in ax.get_xticklabels():
        xtick.set_rotation(90)
    fig = ax.get_figure()
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(os.path.join(output_dir, 'violin.png'))
    fig.savefig(os.path.join(output_dir, 'violin.pdf'))
