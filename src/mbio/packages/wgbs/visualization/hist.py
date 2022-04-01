# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from __future__ import print_function

import logging
import os
import sys

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
    for col_name in df:
        if col_name.endswith('_ratio'):
            series = df[col_name]
            if 'context' not in df:
                series /= 100
            ax = series.hist(bins=100, range=(0, 1))
            ax.set_xlabel('methylation per base')
            ax.set_ylabel('count')
            ax.set_title('Histogram of CpG methylation\nSampleID: {}'.format(col_name[:-6]))
            fig = ax.get_figure()
            fig.tight_layout()
            fig.savefig(os.path.join(output_dir, 'hist.{}.png'.format(col_name[:-6])))
            fig.savefig(os.path.join(output_dir, 'hist.{}.pdf'.format(col_name[:-6])))
            fig.clf()
