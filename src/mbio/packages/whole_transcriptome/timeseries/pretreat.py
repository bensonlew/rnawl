# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'


import os

import numpy as np
import pandas as pd


def main(args):
    params = [args.matrix, args.design, args.level]
    if hasattr(args, 'geneset'):
        params.append(args.geneset)
    print 'INFO: start building data frame'
    df_mat, df_des = table_to_df(*params)
    print 'INFO: start generating delim file for R program'
    delim_mat = os.path.join(args.out_dir, 'matrix.delim.txt')
    delim_des = os.path.join(args.out_dir, 'design.delim.txt')
    df_to_delim(df_mat, delim_mat)
    df_to_delim(df_des, delim_des)
    if os.path.getsize(delim_mat):
        print 'INFO: succeed in exporting file -> {}'.format(delim_mat)
    if os.path.getsize(delim_des):
        print 'INFO: succeed in exporting file -> {}'.format(delim_des)


def table_to_df(matrix_file, design_file, level, geneset_file=None):
    df_mat = pd.read_table(matrix_file)
    df_des = pd.read_table(design_file)
    columns = list(df_des.columns)
    for i in [c for c in columns if 'Unnamed:' in c]:
        columns.remove(i)
    df_des = df_des.reindex(columns, axis=1)
    samples = list(df_des['Sample'])
    df_mat.rename({{'T': 'transcript_id', 'G': 'gene_id'}[level]: 'seq_id'}, axis=1, inplace=True)
    df_mat = df_mat.reindex(['seq_id'] + samples, axis=1)
    if geneset_file:
        seq_ids = [line.strip() for line in open(geneset_file)]
        df_mat = df_mat.query('seq_id in @seq_ids')
    df_mat = df_mat.set_index('seq_id')
    df_mat = np.log10(df_mat + 1)
    df_mat = df_mat.reset_index()
    return df_mat, df_des


def df_to_delim(df, delim):
    lines = ['{}\n'.format('\t'.join(list(df.columns)[1:]))]
    lines.extend('{}\n'.format('\t'.join(str(i) for i in list(record)[1:])) for record in df.to_records())
    open(delim, 'w').writelines(lines)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Pretreat data and design table for MaSigPro.')
    parser.add_argument('-m', dest='matrix', type=str, help='input raw data matrix table', required=True)
    parser.add_argument('-g', dest='geneset', type=str, help='input query geneset list')
    parser.add_argument('-d', dest='design', type=str, help='input design table', required=True)
    parser.add_argument('-l', dest='level', choices=['T', 'G'], help='exp level of matrix', required=True)
    parser.add_argument('-o', dest='out_dir', type=str, help='output result directory', required=True)
    args = parser.parse_args()

    main(args)
