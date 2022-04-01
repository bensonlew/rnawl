# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import pandas as pd
import re
import os

class TableKit:
    def __init__(self, args):
        self.args = args
        self.func = {'combine': self.combine, 'excel': self.excel}[args.subparsers]

    def run(self):
        logging.info('begin of the function ({})'.format(self.func.__name__))
        self.func(self.args)
        logging.info('final of the function ({})'.format(self.func.__name__))

    def combine(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        dfs = [pd.read_csv(line.strip(), sep=sep, index_col=0) for line in open(args.list)]
        for n, df in enumerate(dfs):
            if n == 0:
                index = df.index
                ldf = df
            else:
                ldf = ldf.combine_first(df)
        else:
            if args.method == 'all':
                pass
            elif args.method == 'nona':
                ldf = ldf.dropna()
            elif args.method == 'first':
                ldf = ldf.reindex(index)
            m = re.match(r'(\d+)-', args.fill)
            if m:
                fill = int(m.group(1)) * '-'
            else:
                fill = args.fill
            ldf = ldf.fillna(fill)
            ldf.to_csv(args.output, sep='\t')

    def excel(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        ef = pd.ExcelFile(args.input)
        for sheet_name in ef.sheet_names:
            df = ef.parse(sheet_name, header=None)
            df.to_csv(os.path.join(args.output, '{}.txt'.format(sheet_name)), sep=sep, header=False, index=False)

if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='To be described')
    subparsers = parser.add_subparsers(dest='subparsers', help='sub-command help')

    parser_c = subparsers.add_parser('combine', help='combine several tables')
    parser_c.add_argument('--list', action='store', required=True,
                          help='tables list', dest='list')
    parser_c.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_c.add_argument('--method', action='store', choices=['all', 'nona', 'first'], required=True,
                          help='content style', dest='method')
    parser_c.add_argument('--fill', action='store', required=True,
                          help='default value for null data', dest='fill')
    parser_c.add_argument('--output', action='store', required=True,
                          help='combined table', dest='output')

    parser_e = subparsers.add_parser('excel', help='convert excel to text')
    parser_e.add_argument('--input', action='store', required=True,
                          help='binary excel file', dest='input')
    parser_e.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_e.add_argument('--output', action='store', required=True,
                          help='directory containing results', dest='output')

    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    if len(sys.argv) > 1:
        if sys.argv[1] in ['combine', 'excel']:
            args = parser.parse_args()
        else:
            logging.error('unrecognized command ({})'.format(sys.argv[1]))
            sys.exit(-2)
    else:
        parser.print_help()
        sys.exit(-1)

    inst = TableKit(args)
    inst.run()
