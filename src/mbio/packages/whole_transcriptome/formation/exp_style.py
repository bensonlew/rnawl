# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import os

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
logger = logging.getLogger()


def main(args):
    json_dict = json.load(open(args.json))
    output_dir = args.output

    t_style_df = pd.read_table(json_dict['t_style'], names=['seq_id', 'rna_type', 'rna_kind'],
                               index_col=0, usecols=[0, 2, 3])
    t_style_df['is_new'] = t_style_df['rna_kind'].apply(lambda x: 'true' if x == 'novel' else 'false')
    t_style_df = t_style_df.reindex(['rna_type', 'is_new'], axis=1)

    t_tpm_df = pd.read_table(json_dict['t_tpm'], index_col=0)
    t_tpm_df = t_tpm_df.join(t_style_df)
    t_tpm_df.to_csv(os.path.join(output_dir, os.path.basename(json_dict['t_tpm'])), sep='\t')
    t_fpkm_df = pd.read_table(json_dict['t_fpkm'], index_col=0)
    t_fpkm_df = t_fpkm_df.join(t_style_df)
    t_fpkm_df.to_csv(os.path.join(output_dir, os.path.basename(json_dict['t_fpkm'])), sep='\t')

    g_style_df = pd.read_table(json_dict['g_style'], names=['seq_id', 'rna_type', 'rna_kind'],
                               index_col=0, usecols=[0, 2, 3])
    g_style_df['is_new'] = g_style_df['rna_kind'].apply(lambda x: 'true' if x == 'novel' else 'false')
    g_style_df = g_style_df.reindex(['rna_type', 'is_new'], axis=1)

    g_tpm_df = pd.read_table(json_dict['g_tpm'], index_col=0)
    g_tpm_df = g_tpm_df.join(g_style_df)
    g_tpm_df.to_csv(os.path.join(output_dir, os.path.basename(json_dict['g_tpm'])), sep='\t')
    g_fpkm_df = pd.read_table(json_dict['g_fpkm'], index_col=0)
    g_fpkm_df = g_fpkm_df.join(g_style_df)
    g_fpkm_df.to_csv(os.path.join(output_dir, os.path.basename(json_dict['g_fpkm'])), sep='\t')

    logger.info('succeed in exporting files to {}'.format(output_dir))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate files in old style')
    parser.add_argument('-i', action='store', required=True, help='setting JSON file', dest='json')
    parser.add_argument('-o', action='store', required=True, help='output directory', dest='output')

    args = parser.parse_args()

    main(args)
