# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from __future__ import print_function

import logging
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def func(pct, values):
    return '{:.2f}%\n{}'.format(pct, int(pct / 100 * sum(values)))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: {} <result FILE> <output DIR>'.format(os.path.basename(sys.argv[0])))
    result_fp = sys.argv[1]
    output_dir = sys.argv[2]
    df = pd.read_table(result_fp)
    if 'context' in df:
        for col_name in df:
            if col_name.endswith('_count'):
                value_counts = df[df[col_name] > 5]['context'].value_counts()
                _, ax = plt.subplots(figsize=(6, 6))
                slices, texts, autotexts = ax.pie(value_counts.values,
                                                  autopct=lambda pct: func(pct, value_counts.values),
                                                  textprops={'color': 'w'})
                ax.legend(slices, list(value_counts.index))
                plt.setp(autotexts, size=8, weight='bold')
                ax.set_title('Pie of methylation type\nSampleID: {}'.format(col_name[:-6]))
                fig = ax.get_figure()
                fig.tight_layout()
                fig.savefig(os.path.join(output_dir, 'pie.{}.png'.format(col_name[:-6])))
                fig.savefig(os.path.join(output_dir, 'pie.{}.pdf'.format(col_name[:-6])))
                fig.clf()
