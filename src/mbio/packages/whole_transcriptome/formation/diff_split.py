# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def main(args):
    input_dir = args.input
    relate_file = args.relation
    output_dir = args.output

    relate_df = pd.read_table(relate_file, names=['transcript_id', 'category'], index_col=0, usecols=[0, 2])

    for rna in ('mrna', 'lncrna'):
        if not os.path.isdir(os.path.join(output_dir, rna)):
            os.makedirs(os.path.join(output_dir, rna))

    for fname in os.listdir(input_dir):
        df = pd.read_table(os.path.join(input_dir, fname), index_col=0)
        df = df.join(relate_df)
        df.index.name = 'seq_id'
        mrna_df = df[df['category'] == 'mRNA']
        mrna_df = mrna_df.drop(['category'], axis=1)
        mrna_df.to_csv(os.path.join(output_dir, 'mrna', fname), sep='\t')
        lncrna_df = df[df['category'] == 'lncRNA']
        lncrna_df = lncrna_df.drop(['category'], axis=1)
        lncrna_df.to_csv(os.path.join(output_dir, 'lncrna', fname), sep='\t')
        logger.info('succeed in parsing {}'.format(os.path.join(input_dir, fname)))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Split DE results from longRNA library')
    parser.add_argument('-i', action='store', required=True,
                        help='DE results', metavar='<DIR>', dest='input')
    parser.add_argument('-t', action='store', required=True,
                        help='relation map', metavar='<FILE>', dest='relation')
    parser.add_argument('-o', action='store', required=True,
                        help='output path', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
