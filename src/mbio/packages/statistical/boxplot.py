#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

ap= argparse.ArgumentParser(description='''
boxplot app. need modules: argparse, numpy, pandas, matplotlib, xlwt.
''')
ap.add_argument('data', help='data file, one sample per column, with first row as header and first column as index')
ap.add_argument('-p', '--prefix', help='out prefix', default='test')
args= ap.parse_args()

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot


dat = pd.read_table(args.data, header=0, index_col=0)

d_max = dat.max()
d_min = dat.min()
d_25 = dat.quantile(q=0.25)
d_m = dat.median()
d_75 = dat.quantile(q=0.75)


outlier1 = 3. / 2 * (d_75 - d_25) + d_75
dat_o1 = outlier1 - dat >= 0

outlier2 = d_25 - 3. / 2 * (d_75 - d_25)
dat_o2 = dat - outlier2 >= 0

dat1 = dat[(dat_o1) & (dat_o2)]

outlier_ = dict()

for i in dat.columns:
    outlier = [j for j in dat[~ dat_o1][i] if not np.isnan(j)]
    tmp = [j for j in dat[~ dat_o2][i] if not np.isnan(j)]
    outlier.extend(tmp)
    outlier.sort()
    outlier_[i] = '; '.join([str(j) for j in outlier])

outlier_ = pd.Series(outlier_)

stat = pd.DataFrame([dat1.min(), dat.quantile(q=0.25), dat.median(), dat.quantile(q=0.75), dat1.max(), outlier_])

stat.index = ['min', 'q1', 'median', 'q3', 'max', 'outlier']
stat = stat.T

stat.to_csv('.'.join([args.prefix, 'boxplot.xls']), sep="\t")

cols = dat.columns
dat = dat.values

fig, ax= plot.subplots(figsize=(8, 8))
ax.yaxis.grid(True)
ax.boxplot(dat)

# ax.set_xticks([y+1 for y in range(index.size)])

plot.setp(ax, xticklabels=cols.values)
         
for label in ax.get_xmajorticklabels():
    label.set_rotation(45)

plot.savefig('.'.join([args.prefix, 'boxplot.pdf']))