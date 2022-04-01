#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/12/16 9:59
@file    : LC_measure_drop_qc.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
import pandas as pd
import re
import argparse


# try:
#     _, measure_file, drop_qcs = sys.argv
# except:
#     exit('USAGE: python %s measure_file drop_qc1,drop_qc2' % sys.argv[0] + ' ---drop_qc must be the column of measure_file')


def del_measure_qc(quant_file, del_cols, sort_qc, drop_qc_file):
    measure_file = quant_file
    drop_qcs = del_cols

    with open(measure_file, 'rb') as mr:
        measure_df = pd.read_csv(mr, sep=',', quoting=0, quotechar='"', skiprows=[0, 1], low_memory=False)
    old_columns = measure_df.columns.tolist()

    drop_qcs = drop_qcs.split(',')
    while '' in drop_qcs:
        drop_qcs.remove('')

    # for i in range(len(drop_qcs)):
    #     for c in old_columns:
    #         if c.endswith('_' + drop_qcs[i]):
    #             drop_qcs[i] = c
    #
    # for dq in drop_qcs:
    #     if dq not in measure_df.columns:
    #         exit(dq + ' not in the columns of measure_file')
    yiwai_qc = list()
    for i in range(len(drop_qcs)):
        for c in old_columns:
            if 'QC' in c and '_' in c.split('QC')[-1] and '.' not in c.split('QC')[-1]:
                yiwai_qc.append(c)
            if c.endswith('_' + drop_qcs[i]):
                drop_qcs[i] = c
    drop_qcs += yiwai_qc
    drop_qcs = list(set(drop_qcs))
    for dq in drop_qcs:
        if dq not in measure_df.columns:
            # exit(dq + ' not in the columns of measure_file')
            print(dq + ' not in the columns of measure_file')

    drop_cols = drop_qcs + [qc + '.1' for qc in drop_qcs]
    measure_df = measure_df.drop(columns=drop_cols)
    re_qc = re.compile('QC\d')
    qc_cols = [col for col in measure_df.columns if re_qc.search(col)]

    with open(measure_file) as mr, open(drop_qc_file, 'w') as mw:
        head1 = mr.readline().strip().split(',')
        head2 = mr.readline().strip().split(',')
        relas = list(zip(old_columns, head1, head2))
        head1 = list()
        head2 = list()
        head3 = list()
        for n, (a, b, c) in enumerate(relas):
            a = a.strip('"')
            if a in drop_cols:
                if b or c:
                    relas[n + 1] = (relas[n + 1][0], b, c)
                continue
            head1.append(b)
            head2.append(c)
            if a.endswith('.1'):
                a = a.split('.1')[0]
            if sort_qc:
                if a in qc_cols:
                    a = a.split('QC')[0] + 'QC' + str(qc_cols.index(a) + 1).zfill(2)
            head3.append('"' + a + '"')
        mw.write('\n'.join([','.join(head1), ','.join(head2), ','.join(head3)]) + '\n')

    measure_df.to_csv(drop_qc_file, sep=',', quoting=0, quotechar='"', header=False, index=False, mode='a')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-quant", "--quant", dest="quant", required=True, type=str, help="")
    parser.add_argument("-del_qcs", "--del_qcs",dest="del_qcs", required=True, type=str, help="")
    parser.add_argument("-sort_qc", "--sort_qc",dest="sort_qc", required=True, type=int, help="")
    parser.add_argument("-drop_qc", "--drop_qc",dest="drop_qc", required=True, type=str, help="")
    args = parser.parse_args()
    del_measure_qc(args.quant, args.del_qcs, args.sort_qc, args.drop_qc)
