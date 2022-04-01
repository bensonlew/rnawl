# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from __future__ import print_function

import logging
import os
import sys

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def parse_gtf(filepath):
    m_dict = {'gene': 'Gene body', 'exon': 'Exon', 'CDS': 'CDS'}
    lines = open(filepath).readlines()
    anno_dict = {i: set() for i in range(1, int(lines[-1].split('\t')[4]) + 1)}
    for line in lines:
        eles = line.strip().split('\t')
        for i in range(int(eles[3]), int(eles[4]) + 1):
            if eles[2] in m_dict:
                anno_dict[i].add(m_dict[eles[2]])
        if eles[2] == 'gene':
            __x = int(eles[3]) if int(eles[3]) - 600 > 0 else 1
            __y = int(eles[4]) if int(eles[4]) - 600 > 0 else 1
            for j in range(__x, __y):
                anno_dict[j].add('Promoter')
    for i in anno_dict:
        if len(anno_dict[i]) == 1 and 'Promoter' in anno_dict[i]:
            anno_dict[i].add('Intergenic region')
    return anno_dict


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('usage: {} <gtf FILE> <result FILE> <output DIR>'.format(os.path.basename(sys.argv[0])))
    gtf_fp = sys.argv[1]
    result_fp = sys.argv[2]
    output_dir = sys.argv[3]
    logger.debug('Parse {}'.format(gtf_fp))
    a_dict = parse_gtf(gtf_fp)
    src_df = pd.read_table(result_fp)
    data_dict = {}
    logger.debug('Init data dictionary')
    for col_name in src_df.columns:
        if col_name.endswith('_count'):
            data_dict[col_name[:-6]] = {
                k0: {k1: 0 for k1 in ('Promoter', 'Gene body', 'Intergenic region', 'Exon', 'CDS', 'Intron', 'Total')}
                for k0 in ('CG', 'CHG', 'CHH')}
    logger.debug('Start iteration')
    for _, row in src_df.iterrows():
        context = row['context'] if 'context' in row else 'CG'
        for idx in row.index:
            if idx.endswith('_count') and row[idx] > 5:
                sp_name = idx[:-6]
                if row['pos'] in a_dict:
                    for feature in a_dict[row['pos']]:
                        data_dict[sp_name][context][feature] += 1
                else:
                    data_dict[sp_name][context]['Intergenic region'] += 1
                data_dict[sp_name][context]['Total'] += 1
    for sp_name, data in data_dict.items():
        df = pd.DataFrame(data)
        df.index.name = sp_name
        if 'context' not in src_df:
            df = df.reindex(['CG'], axis=1)
        df.to_csv(os.path.join(output_dir, 'bar.{}.txt'.format(sp_name)), sep='\t')
        ax = df.plot(kind='bar')
        if 'context' not in src_df:
            ax.legend_.remove()
        fig = ax.get_figure()
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, 'bar.{}.png'.format(sp_name)))
        fig.savefig(os.path.join(output_dir, 'bar.{}.pdf'.format(sp_name)))
