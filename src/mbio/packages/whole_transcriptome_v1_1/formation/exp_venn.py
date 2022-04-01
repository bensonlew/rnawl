# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

import logging

import numpy as np
import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    exp_matrix, group_table, threshold, venn_table = args.exp, args.group, args.threshold, args.venn
    exp_df, group_dict = parse_input(exp_matrix, group_table)
    group_exp_df = process_exp_df(exp_df, group_dict, log10=False)
    if export_venn_data(group_exp_df, threshold, venn_table):
        logging.info('succeed in exporting venn data to {}'.format(venn_table))


def export_venn_data(exp_df, threshold, venn_table):
    lines = ['name\tseqs\n']
    for name in exp_df.columns:
        seqs = ','.join(set(exp_df.index[exp_df[name] >= threshold]))
        lines.append('{}\t{}\n'.format(name, seqs))
    else:
        open(venn_table, 'w').writelines(lines)
        return True


def parse_input(exp_matrix, group_table):
    exp_df = pd.read_table(exp_matrix, index_col=0)
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
    else:
        exp_df = exp_df.reindex(sample_list, axis=1)
    return exp_df, group_dict


def process_exp_df(exp_df, group_dict=None, log10=True):
    if group_dict:
        group_exp_dfs = list()
        for group in group_dict:
            group_exp_df = exp_df.reindex(group_dict[group], axis=1).mean(axis=1)
            group_exp_df.name = group
            group_exp_dfs.append(group_exp_df)
        else:
            exp_df = pd.concat(group_exp_dfs, axis=1)
    if log10:
        exp_df = np.log10(exp_df)
    return exp_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate venn table')
    parser.add_argument('-i', action='store', required=True,
                        help='expression matrix', metavar='<FILE>', dest='exp')
    parser.add_argument('-g', action='store', required=True,
                        help='group table', metavar='<FILE>', dest='group')
    parser.add_argument('-t', action='store', type=float, required=True,
                        help='group expression threshold', metavar='<FLOAT>', dest='threshold')
    parser.add_argument('-o', action='store', required=True,
                        help='venn table', metavar='<FILE>', dest='venn')

    args = parser.parse_args()

    main(args)
