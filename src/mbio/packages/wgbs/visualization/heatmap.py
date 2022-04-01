# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from __future__ import print_function

import logging
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('usage: {} <dml.result FILE> <result FILE> <output DIR>'.format(os.path.basename(sys.argv[0])))
    dml_result_fp = sys.argv[1]
    result_fp = sys.argv[2]
    output_dir = sys.argv[3]
    dml_result_df = pd.read_table(dml_result_fp)
    if dml_result_df.shape[0]:
        logger.debug('Found DML, take DML data')
        df = dml_result_df
        fig_name = 'heatmap.dml'
    else:
        logger.debug('Found no DML, take all data')
        df = pd.read_table(result_fp)
        fig_name = 'heatmap'
    logger.debug('Prepare data')
    df['idx'] = df['chrom'] + df['pos'].astype(str) + df['strand']
    df.set_index('idx', inplace=True)
    df = df.reindex(sorted(col_name for col_name in df.columns if col_name.endswith('_ratio')), axis=1)
    df.rename({col_name: col_name.split('-')[2] for col_name in df.columns}, axis=1, inplace=True)
    df.index.name = None
    logger.debug('Draw heatmap')
    _, ax = plt.subplots(figsize=(6, 18))
    sns.set(style='white')
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.heatmap(df, cmap=cmap, ax=ax)
    for xtick in ax.xaxis.get_major_ticks():
        xtick.tick1line.set_markersize(0)
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, '{}.png'.format(fig_name)))
    fig.savefig(os.path.join(output_dir, '{}.pdf'.format(fig_name)))
