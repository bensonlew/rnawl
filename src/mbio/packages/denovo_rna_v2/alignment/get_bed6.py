# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging

import pandas as pd

from mbio.packages.denovo_rna_v2.functions import watcher


@watcher
def main(args):
    idf = pd.read_table(args.input,
                        names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                               'sstart', 'send', 'evalue', 'bitscore'])
    parse = lambda i, row: (
        row['sseqid'], min(row['sstart'] - 1, row['send']), max(row['sstart'] - 1, row['send']), row['qseqid'], i,
        '+' if row['sstart'] < row['send'] else '-')
    data = [parse(i, row) for i, row in idf.iterrows() if row['evalue'] <= args.evalue]
    odf = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    odf.to_csv(args.output, sep='\t', header=False, index=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate BED6 from BLAST M8 tabular')

    parser.add_argument('--input', action='store', required=True,
                        help='blast m8 tabular', metavar='<FILE>', dest='input')
    parser.add_argument('--evalue', action='store', type=float, required=True,
                        help='threshold', metavar='<STR>', dest='evalue')
    parser.add_argument('--output', action='store', required=True,
                        help='6 columns bed', metavar='<FILE>', dest='output')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    main(args)
