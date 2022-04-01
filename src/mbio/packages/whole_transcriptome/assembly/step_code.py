# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import pickle

import pandas as pd
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    fasta = args.fasta
    tmap = args.tmap
    output_dir = args.output

    seq_len_list = [len(record) for record in SeqIO.parse(fasta, 'fasta')]
    step_data = list()
    for step in [200, 300, 600, 1000]:
        step_dict = count_step(step, seq_len_list)
        step_data.append({'step': step, 'step_data': step_dict})
    else:
        step_file = os.path.join(output_dir, 'step.pk')
        pickle.dump(step_data, open(step_file, 'w'))
        logging.info('succeed in exporting {}'.format(step_file))
        code_data = class_code(fasta, tmap)
        code_file = os.path.join(output_dir, 'code.pk'.format(step))
        pickle.dump(code_data, open(code_file, 'w'))
        logging.info('succeed in exporting {}'.format(code_file))


def count_step(step, seq_len_list, y=10):
    step_data = list()
    step_dict = dict((step * i, int()) for i in range(1, y))
    over_step = 0
    for seq_len in seq_len_list:
        for upper_limit in sorted(step_dict.keys()):
            if seq_len <= upper_limit:
                step_dict[upper_limit] += 1
                break
        else:
            over_step += 1
    else:
        for upper_limit in sorted(step_dict.keys()):
            step_data.append({'{}~{}'.format(upper_limit - step + 1, upper_limit): step_dict[upper_limit]})
        else:
            step_data.append({'>{}'.format(step * (y - 1)): over_step})
            step_data.append({'total': len(seq_len_list)})
    return step_data


def class_code(fasta, tmap):
    # cc_series = pd.read_table(tmap, usecols=['class_code'])['class_code']
    # vc_series = cc_series.value_counts()
    # vc_df = vc_series.to_frame()
    # cn_df = vc_df.reset_index().rename({'index': 'class_code', 'class_code': 'num'}, axis=1)
    # code_data = cn_df.to_dict('r')
    df = pd.read_table(tmap)
    df = df.reindex(['qry_id', 'class_code'], axis=1)
    df = df.set_index('qry_id')
    cc2nu_dict = {cc: 0 for cc in df['class_code'].unique()}
    id2cc_dict = df.to_dict()['class_code']
    for record in SeqIO.parse(fasta, 'fasta'):
        cc = id2cc_dict.get(record.id, '=')
        cc2nu_dict[cc] += 1
    code_data = [{'class_code': class_code, 'num': num} for class_code, num in cc2nu_dict.items()]
    return code_data


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate data for dumping into mongo')
    parser.add_argument('-f', action='store', required=True,
                        help='all transcripts FASTA', metavar='<FILE>', dest='fasta')
    parser.add_argument('-t', action='store', required=True,
                        help='compare TMAP', metavar='<FILE>', dest='tmap')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
