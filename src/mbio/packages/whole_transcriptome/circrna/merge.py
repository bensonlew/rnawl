# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import pandas as pd


def main(args):
    f1 = pd.read_table(args.find_circ)
    f1['chr'] = f1['chr'].astype(str)
    f2 = pd.read_table(args.CIRI2)
    f2['chr'] = f2['chr'].astype(str)
    f3 = pd.merge(f1, f2, how='outer', on=['circRNA_name', 'chr', 'circRNA_start', 'circRNA_end', 'strand'],
                  indicator=True)
    f3['_merge'].replace('left_only', 'find_circ', inplace=True)
    f3['_merge'].replace('right_only', 'CIRI2', inplace=True)
    f3.sort_values(by=['chr', 'circRNA_start', 'circRNA_end'], ascending=True, inplace=True)
    f3.to_csv(args.merge, index=False, sep='\t')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='merge find_circ&CIRC2')
    parser.add_argument('-f', action='store', required=True,
                        help='find_circ file', dest='find_circ')
    parser.add_argument('-c', action='store', required=True,
                        help='CIRC2 file', dest='CIRI2')
    parser.add_argument('-o', action='store', required=True,
                        help='merge ', dest='merge')

    args = parser.parse_args()
    main(args)
