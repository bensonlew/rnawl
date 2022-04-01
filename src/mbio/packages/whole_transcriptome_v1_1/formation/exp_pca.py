# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

import logging
import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    exp_matrix, group_table, mean, output_dir = args.exp, args.group, args.mean, args.output
    exp_df, group_dict = parse_input(exp_matrix, group_table)
    exp_df = process_exp_df(exp_df, group_dict if mean else None)
    if export_pca_data(exp_df, output_dir):
        logging.info('succeed in exporting corr data to {}'.format(output_dir))


def export_pca_data(exp_df, output_dir):
    X = exp_df.T
    pca = PCA()
    pca.fit(X)
    evr_lines = list()
    cumsum_ratio = 0.0
    n_components = 2
    for i, ratio in enumerate(pca.explained_variance_ratio_):
        evr_lines.append('PC{}\t{}\n'.format(i + 1, ratio))
        cumsum_ratio += ratio
        if cumsum_ratio >= 0.95:
            n_components = i + 1
            break
    X_transformed = pca.transform(X)
    pca_df = pd.DataFrame(X_transformed, index=pd.Index(X.index, name='sample'))
    pca_df = pca_df.iloc[:, :n_components]
    pca_df.columns = ['PC{}'.format(n + 1) for n in range(n_components)]
    pca_df = pca_df.round(2)
    pca_df.to_csv(os.path.join(output_dir, 'pca.txt'), sep='\t')
    open(os.path.join(output_dir, 'explained_variance_ratio.txt'), 'w').writelines(evr_lines)
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
        exp_df = np.log10(exp_df + 1)
    return exp_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate pca table')
    parser.add_argument('-i', action='store', required=True,
                        help='expression matrix', metavar='<FILE>', dest='exp')
    parser.add_argument('-g', action='store', required=True,
                        help='group table', metavar='<FILE>', dest='group')
    parser.add_argument('-m', '--mean', action='store_true',
                        help='take group mean of values', dest='mean')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
