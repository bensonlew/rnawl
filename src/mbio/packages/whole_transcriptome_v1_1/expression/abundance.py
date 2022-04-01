# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import os
import sys

import numpy as np
import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def main(args):
    json_dict = json.load(open(args.json))
    logger.info('succeed in parsing arguments as follow ->\n{}'.format(json_dict))

    m_table = json_dict['m_table']
    i_table = json_dict['i_table']
    t_table = json_dict['t_table']
    g_table = json_dict['g_table']

    m_df = pd.read_table(m_table, names=['transcript_id', 'gene_id'], index_col=0)
    i_df = pd.read_table(i_table, index_col=0)
    i_df.index.name = 'transcript_id'
    j_df = i_df.join(m_df).reset_index().set_index(['gene_id', 'transcript_id'])
    denom = j_df['eff_length'].sum()
    j_df['fpkm'] = j_df['est_counts'] / j_df['eff_length'] / denom * 1e9
    j_df = j_df.rename({'est_counts': 'count'}, axis=1).reindex(['tpm', 'fpkm', 'count'], axis=1)
    t_df = j_df.reset_index().set_index('transcript_id').drop('gene_id', axis=1)
    g_df = j_df.sum(level='gene_id')
    t_df.to_csv(t_table, sep='\t')
    g_df.to_csv(g_table, sep='\t')

    for name in ('t_table', 'g_table'):
        if os.path.isfile(json_dict[name]):
            logger.info('succeed in exporting {}'.format(json_dict[name]))


def get_indices(target_ids, transcript_ids):
    arr = np.array(target_ids)
    find = lambda x: np.where(arr == x)[0][0]
    indices = map(find, transcript_ids)
    assert np.array(target_ids) == np.array(transcript_ids)[indices]


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Build abundance tables')
    parser.add_argument('json', action='store', help='setting file in JSON format')

    args = parser.parse_args()

    if hasattr(args, 'json'):
        main(args)
    else:
        parser.print_help()
        sys.exit(-1)
