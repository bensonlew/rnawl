# -*- coding:utf-8 -*-

import glob
import logging
import os
import shutil
import sys

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

if __name__ == '__main__':
    dirs_path = sys.argv[1]
    if os.path.isdir(os.path.join(dirs_path, 'output')):
        shutil.rmtree(os.path.join(dirs_path, 'output'))
    dir_paths = glob.glob(os.path.join(dirs_path, '*'))
    df_dict = dict()
    max_sp_name, max_n_rows = str(), int()
    for path in dir_paths:
        if not os.path.isdir(path):
            continue
        sp_name = os.path.basename(path)
        txt_fp = os.path.join(os.path.join(path, '{}.txt'.format(sp_name)))
        logger.debug('Parse {}'.format(txt_fp))
        txt_df = pd.read_table(txt_fp)
        txt_df['idx'] = txt_df['chr'] + '_' + txt_df['pos'].astype(str) + '_' + txt_df['strand'] + '_' + txt_df[
            'context']
        txt_df.rename({'chr': 'chrom', 'ratio': '{}_ratio'.format(sp_name), 'C_count': '{}_count'.format(sp_name)},
                      axis=1, inplace=True)
        txt_df = txt_df.reindex(['idx', '{}_ratio'.format(sp_name), '{}_count'.format(sp_name)], axis=1).set_index(
            'idx')
        logger.debug('Show data frame shape of {}: {}'.format(sp_name, txt_df.shape))
        df_dict[sp_name] = txt_df
        if df_dict[sp_name].shape[0] > max_n_rows:
            max_sp_name = sp_name
            max_n_rows = df_dict[sp_name].shape[0]
    logger.debug('Adjust index of data frame')
    for sp_name, df in df_dict.items():
        df_dict[sp_name] = df.reindex(df_dict[max_sp_name].index)
        logger.debug('Show data frame shape of {}: {}'.format(sp_name, df_dict[sp_name].shape))
    logger.debug('Merge and reshape combined data frame')
    df = pd.concat(df_dict.values(), axis=1)
    for col_name in df:
        if col_name.endswith('_count'):
            df[col_name] = df[col_name].fillna(0).astype(int)
    columns = list(df.columns)
    chrom_list, pos_list, strand_list, context_list = list(), list(), list(), list()
    for eles in df.index.str.split('_'):
        chrom_list.append(eles[0])
        pos_list.append(eles[1])
        strand_list.append(eles[2])
        context_list.append(eles[3])
    df['chrom'] = chrom_list
    df['pos'] = pos_list
    df['strand'] = strand_list
    df['context'] = context_list
    columns = ['chrom', 'pos', 'strand', 'context'] + columns
    df = df.reindex(columns, axis=1)
    output_dir = os.path.join(dirs_path, 'output')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    logger.debug('Export {}'.format(os.path.join(output_dir, 'result.txt')))
    df.to_csv(os.path.join(output_dir, 'result.txt'), sep='\t', index=False)
