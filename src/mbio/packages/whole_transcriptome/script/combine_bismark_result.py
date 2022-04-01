# -*- coding:utf-8 -*-

import glob
import gzip
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
    for path in dir_paths:
        if not os.path.isdir(path):
            continue
        sp_name = os.path.basename(path)
        cov_fp = os.path.join(os.path.join(path, 'result/{}_pe.deduplicated.bismark.cov'.format(sp_name)))
        if not os.path.isfile(cov_fp) and os.path.isfile('{}.gz'.format(cov_fp)):
            logger.debug('Unzip {}.gz'.format(cov_fp))
            with gzip.GzipFile('{}.gz'.format(cov_fp)) as gz_file:
                open(cov_fp, 'wb+').write(gz_file.read())
        logger.debug('Parse {}'.format(cov_fp))
        cov_df = pd.read_table(cov_fp, header=None)
        cov_df[6] = cov_df[0] + cov_df[1].astype(str)
        cov_df.rename({0: 'chrom', 1: 'pos', 3: 'ratio', 4: 'count', 6: 'idx'}, axis=1, inplace=True)
        cov_df = cov_df.reindex(['idx', 'chrom', 'pos', 'ratio', 'count'], axis=1).set_index('idx')
        txt_fp = os.path.join(os.path.join(path, 'result/{}_pe.deduplicated.CpG_report.txt'.format(sp_name)))
        logger.debug('Parse {}'.format(txt_fp))
        txt_df = pd.read_table(txt_fp, header=None)
        txt_df[7] = txt_df[0] + txt_df[1].astype(str)
        txt_df.rename({2: 'strand', 7: 'idx'}, axis=1, inplace=True)
        txt_df = txt_df.reindex(['idx', 'strand'], axis=1).set_index('idx')
        logger.debug('Generate merged data frame for {}'.format(sp_name))
        df = cov_df.join(txt_df)
        df = df[pd.notna(df['strand'])]
        output_dir = os.path.join(path, 'output')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        df.to_csv(os.path.join(output_dir, '{}.detail.txt'.format(sp_name)), sep='\t', index=False)
        df_dict[sp_name] = df
    max_sp_name, max_n_rows = str(), int()
    for sp_name, df in df_dict.items():
        df['idx'] = df['chrom'] + '_' + df['pos'].astype(str) + '_' + df['strand']
        df.rename({'ratio': '{}_ratio'.format(sp_name), 'count': '{}_count'.format(sp_name)}, axis=1, inplace=True)
        df = df.reindex(['idx', '{}_ratio'.format(sp_name), '{}_count'.format(sp_name)], axis=1).drop_duplicates()
        df_dict[sp_name] = df.set_index('idx')
        logger.debug('Show data frame shape of {}: {}'.format(sp_name, df_dict[sp_name].shape))
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
    chrom_list, pos_list, strand_list = list(), list(), list()
    for eles in df.index.str.split('_'):
        chrom_list.append(eles[0])
        pos_list.append(eles[1])
        strand_list.append(eles[2])
    df['chrom'] = chrom_list
    df['pos'] = pos_list
    df['strand'] = strand_list
    columns = ['chrom', 'pos', 'strand'] + columns
    df = df.reindex(columns, axis=1)
    output_dir = os.path.join(dirs_path, 'output')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    logger.debug('Export {}'.format(os.path.join(output_dir, 'result.txt')))
    df.to_csv(os.path.join(output_dir, 'result.txt'), sep='\t', index=False)
