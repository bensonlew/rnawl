# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

import logging
import os

import numpy as np
import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    exp_matrix, group_table, method, mean, log10, output_dir = args.exp, args.group, args.method, args.mean, args.log, args.output
    exp_df, group_dict = parse_input(exp_matrix, group_table)
    exp_df = process_exp_df(exp_df, group_dict if mean else None, log10)
    if export_corr_data(exp_df, method, output_dir):
        logging.info('succeed in exporting corr data to {}'.format(output_dir))


def export_corr_data(exp_df, method, output_dir):
    corr_df = exp_df.corr(method=method)
    corr_df.index.name = 'sample'
    corr_df.to_csv(os.path.join(output_dir, 'corr.txt'), sep='\t')
    exp_df.to_csv(os.path.join(output_dir, 'delim.txt'), sep='\t')
    return True


def parse_input(exp_matrix, group_table):
    exp_df = pd.read_table(exp_matrix, index_col=0)
    index_name = exp_df.index.name
    exp_df = exp_df.reset_index().drop_duplicates().set_index(index_name)
    group_dict = dict()
    sample_list = list()
    for line in open(group_table):
        if line.strip() and line[0] != '#':
            sample, group = line.strip().split('\t')
            sample_list.append(sample)
            if group in group_dict:
                group_dict[group].append(sample)
            else:
                group_dict[group] = [sample]
    exp_df = exp_df.reindex(sample_list, axis=1)
    return exp_df, group_dict


def process_exp_df(exp_df, group_dict=None, log10=False):
    if group_dict:
        group_exp_dfs = list()
        for group in group_dict:
            group_exp_df = exp_df.reindex(group_dict[group], axis=1).mean(axis=1)
            group_exp_df.name = group
            group_exp_dfs.append(group_exp_df)
        else:
            exp_df = pd.concat(group_exp_dfs, axis=1)
    if log10:
        exp_df = np.log10(exp_df + 1)
    return exp_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate corr table')
    parser.add_argument('-i', action='store', required=True,
                        help='expression matrix', metavar='<FILE>', dest='exp')
    parser.add_argument('-g', action='store', required=True,
                        help='group table', metavar='<FILE>', dest='group')
    parser.add_argument('-m', action='store', choices=['pearson', 'kendall', 'spearman'], required=True,
                        help='pairwise correlation method', dest='method')
    parser.add_argument('-t', '--mean', action='store_true',
                        help='take group mean of values', dest='mean')
    parser.add_argument('-l', '--log', action='store_true',
                        help='take log 10 of values', dest='log')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
