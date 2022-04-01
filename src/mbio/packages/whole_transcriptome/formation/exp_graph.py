# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

import logging
import os
import pickle

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    exp_matrix, group_table, output_dir, log10 = args.exp, args.group, args.output, args.log
    exp_df, group_dict = parse_input(exp_matrix, group_table)
    if export_box_data(exp_df, group_dict, output_dir, log10):
        logging.info('succeed in exporting box data to {}'.format(output_dir))
    if export_density_data(exp_df, group_dict, output_dir, log10):
        logging.info('succeed in exporting density data to {}'.format(output_dir))
    if export_volin_data(exp_df, group_dict, output_dir, log10):
        logging.info('succeed in exporting volin data to {}'.format(output_dir))


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


def export_box_data(exp_df, group_dict, output_dir, log10):
    sample_exp_df = process_exp_df(exp_df, log10=log10)
    sample_exp_df = sample_exp_df[sample_exp_df.sum(axis=1) > 1e-3]
    if len(sample_exp_df) > 1:
        sample_box_data = get_box_data(sample_exp_df)
    else:
        sample_box_data = dict()
    pickle.dump(sample_box_data, open(os.path.join(output_dir, 'sample_box_data.pk'), 'w'))
    group_exp_df = process_exp_df(exp_df, group_dict, log10=log10)
    group_exp_df = group_exp_df[group_exp_df.sum(axis=1) > 1e-3]
    if len(group_exp_df) > 1:
        group_box_data = get_box_data(group_exp_df)
    else:
        group_box_data = dict()
    pickle.dump(group_box_data, open(os.path.join(output_dir, 'group_box_data.pk'), 'w'))
    return True


def get_box_data(exp_df):
    box_data = list()
    for column in exp_df.columns:
        exp_series = exp_df[column]
        exp_series = exp_series[exp_series != 0]
        des_series = exp_series.describe()
        doc_dict = {
            'sample': column,
            'max': des_series['max'],
            'q3': des_series['75%'],
            'median': des_series['50%'],
            'q1': des_series['25%'],
            'min': des_series['min']
        }
        box_data.append(doc_dict)
    else:
        return box_data


def export_density_data(exp_df, group_dict, output_dir, log10):
    sample_exp_df = process_exp_df(exp_df, log10=log10)
    sample_exp_df = sample_exp_df[sample_exp_df.sum(axis=1) > 1e-3]
    if len(sample_exp_df) > 1:
        sample_density_data = get_density_data(sample_exp_df)
    else:
        sample_density_data = dict()
    pickle.dump(sample_density_data, open(os.path.join(output_dir, 'sample_density_data.pk'), 'w'))
    group_exp_df = process_exp_df(exp_df, group_dict, log10=log10)
    group_exp_df = group_exp_df[group_exp_df.sum(axis=1) > 1e-3]
    if len(group_exp_df) > 1:
        group_density_data = get_density_data(group_exp_df)
    else:
        group_density_data = dict()
    pickle.dump(group_density_data, open(os.path.join(output_dir, 'group_density_data.pk'), 'w'))
    return True


def get_density_data(exp_df):
    density_data = list()
    for column in exp_df.columns:
        exp_series = exp_df[column]
        exp_series = exp_series[exp_series != 0]
        kernel = stats.gaussian_kde(exp_series)
        points = np.linspace(exp_series.min(), exp_series.max(), num=1000, endpoint=False)
        result = kernel(points)
        data = pd.DataFrame({'log2exp': points, 'density': result}).to_dict('r')
        doc_dict = {
            'sample': column,
            'data': data
        }
        density_data.append(doc_dict)
    else:
        return density_data


def export_volin_data(exp_df, group_dict, output_dir, log10):
    sample_exp_df = process_exp_df(exp_df, log10=log10)
    sample_exp_df = sample_exp_df[sample_exp_df.sum(axis=1) > 1e-3]
    if len(sample_exp_df) > 1:
        sample_volin_data = get_volin_data(sample_exp_df)
    else:
        sample_volin_data = dict()
    pickle.dump(sample_volin_data, open(os.path.join(output_dir, 'sample_volin_data.pk'), 'w'))
    group_exp_df = process_exp_df(exp_df, group_dict, log10=log10)
    group_exp_df = group_exp_df[group_exp_df.sum(axis=1) > 1e-3]
    if len(group_exp_df) > 1:
        group_volin_data = get_volin_data(group_exp_df)
    else:
        group_volin_data = dict()
    pickle.dump(group_volin_data, open(os.path.join(output_dir, 'group_volin_data.pk'), 'w'))
    return True


def get_volin_data(exp_df):
    frac_exp_df = exp_df.sample(frac=0.7).reset_index()
    volin_data = frac_exp_df.to_dict('r')
    return volin_data


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

    parser = argparse.ArgumentParser(description='Generate files for dumping data into mongo')
    parser.add_argument('-i', action='store', required=True,
                        help='expression matrix', metavar='<FILE>', dest='exp')
    parser.add_argument('-g', action='store', required=True,
                        help='group table', metavar='<FILE>', dest='group')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')
    parser.add_argument('-l', action='store_true', help='take log', dest='log')

    args = parser.parse_args()

    main(args)
