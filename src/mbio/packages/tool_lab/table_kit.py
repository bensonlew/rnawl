# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os
import re
import subprocess
import numpy as np
import pandas as pd
from sklearn import preprocessing

class TableKit(object):
    def __init__(self, args):
        self.args = args
        self.func = {
            'combine': self.combine,
            'excel': self.excel,
            'transpose': self.transpose,
            'fillna': self.fillna,
            'standard': self.standard,
            'filter': self.filter
        }[args.subparsers]

    def run(self):
        logging.info('begin of the function ({})'.format(self.func.__name__))
        self.func(self.args)
        logging.info('final of the function ({})'.format(self.func.__name__))

    def combine(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        table_list = args.list.split(',')
        dfs = [pd.read_csv(line.strip(), sep=sep, index_col=0) for line in table_list]
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
            if args.fill == 'doubleline':
                fill = '--'
            ldf = ldf.fillna(fill)
            ldf.to_csv(args.output, sep='\t')

    def excel(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': ' '}[args.sep]
        ef = pd.ExcelFile(args.input)
        for sheet_name in ef.sheet_names:
            df = ef.parse(sheet_name, header=None)
            df.to_csv(os.path.join(args.output, '{}.txt'.format(sheet_name)), sep=sep, header=False, index=False)

    def transpose(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        df = pd.read_csv(args.input, sep=sep,index_col=0)
        tdf = df.T
        tdf.to_csv(args.output, sep='\t', header=True,index=True)

    def fillna(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        na_values = [''] if args.nav == 'blank' else [args.nav]
        df = pd.read_csv(args.input, sep=sep, na_values=na_values, keep_default_na=False, index_col=0)
        if args.rmp:
            thresh = df.shape[1] * (args.rmp / 100)
            df.dropna(thresh=thresh, inplace=True)
        else:
            pass
        if args.method in ['missForest', 'kNN', 'PPCA', 'BPCA', 'nipals']:
            r_input = os.path.join(os.path.dirname(args.output), 'for.r.tmp')
            df.to_csv(r_input, sep='\t', index=True)
            cmd = '{} {} -i {} -m {} -o {}'.format(
                args.r_interpreter, args.r_script, r_input, args.method.lower(), args.output)
            proc = subprocess.Popen(
                cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            retcode = proc.wait()
            outs, errs = proc.communicate()
            if retcode:
                msg = '\n'.join(
                    ('fail to excecute {}'.format(cmd), 'STDOUT: {}'.format(outs), 'STDERR: {}'.format(errs)))
                raise Exception(msg)
            else:
                logging.info('succeed in excecuting command ({})'.format(cmd))
                logging.debug(outs)
                logging.error(errs)
        else:
            if args.method == 'min':
                value = min(df.min())
            elif args.method == 'otMin':
                value = min(df.min()) / 10.0
            elif args.method == 'halfMin':
                value = min(df.min()) / 2.0
            elif args.method == 'colMin':
                value = {c: df[c].min() for c in df.columns}
            elif args.method == 'colMean':
                value = {c: df[c].mean() for c in df.columns if df[c].dtype != object}
            elif args.method == 'colMedian':
                value = {c: df[c].median() for c in df.columns if df[c].dtype != object}
            df = df.fillna(value)
            df.to_csv(args.output, sep='\t', index=True)

    def standard(self, args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        df = pd.read_csv(args.input, sep=sep, index_col=0)
        if args.observe == 'none':
            pass
        elif args.observe == 'sum':
            df /= df.sum()
        elif args.observe == 'median':
            df /= df.median()
        elif args.observe == 'mean':
            df /= df.mean()
        elif args.observe == 'refrow':
            refrow = df[args.guard].copy()
            df = pd.concat([df[i] / refrow.rename(df[i].name) for i in df.columns], axis=1)
        elif args.observe == 'refcol':
            refcol = df.loc[args.guard].copy()
            df /= refcol
        if args.feature == 'none':
            pass
        elif args.feature == 'standard_scale':
            df = pd.DataFrame(preprocessing.scale(df, with_std=False, axis=1), index=df.index, columns=df.columns)
        elif args.feature == 'minmax_scale':
            df = pd.DataFrame(preprocessing.minmax_scale(df, axis=1), index=df.index, columns=df.columns)
        elif args.feature == 'maxabs_scale':
            # df = pd.DataFrame(preprocessing.maxabs_scale(df, axis=1), index=df.index, columns=df.columns)
            df_T = (df.T - df.T.mean())/np.sqrt(np.std(df.T))
            df = df_T.T
        elif args.feature == 'normalize':
            df = pd.DataFrame(preprocessing.scale(df, axis=1), index=df.index, columns=df.columns)
        df[df==-0]=0
        df.to_csv(args.output, sep='\t')

    def filter(self,args):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[args.sep]
        if args.table_head is True:
            df = pd.read_csv(args.input, sep=sep, index_col=False, header=0)
        else:
            df = pd.read_csv(args.input, sep=sep, index_col=False, header=None)
        filter_list = list()
        with open(args.info, 'r') as filtered:
            for line in filtered.readlines():
                filter_list.append(line.strip())
        if args.lower is True:
            filter_list = [x.lower() for x in filter_list]
        else:
            pass
        if args.contain is False:
            if args.lower is True:
                df_filter = df[df.iloc[:, int(args.number)-1].str.lower().isin(filter_list)]
            else:
                df_filter = df[df.iloc[:, int(args.number) - 1].isin(filter_list)]
        else:
            if args.lower is True:
                pattern = '|'.join(filter_list)
                df_filter = df[df.iloc[:, int(args.number)-1].str.contains(pattern, case=False)]
            else:
                pattern = '|'.join(filter_list)
                df_filter = df[df.iloc[:, int(args.number)-1].str.contains(pattern, case=True)]
        print df_filter
        df_filter.to_csv(args.output, sep='\t', index=False)




if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Tool kit for manipulating text table file')
    subparsers = parser.add_subparsers(dest='subparsers', help='sub-command help')

    parser_c = subparsers.add_parser('combine', help='combine several tables')
    parser_c.add_argument('--list', action='store', required=True,
                          help='tables list', dest='list')
    parser_c.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_c.add_argument('--method', action='store', choices=['all', 'nona', 'first'], required=True,
                          help='content style', dest='method')
    parser_c.add_argument('--fill', action='store', required=True,
                          help='default value for NA data', dest='fill')
    parser_c.add_argument('--output', action='store', required=True,
                          help='combined table', dest='output')

    parser_e = subparsers.add_parser('excel', help='convert excel to text')
    parser_e.add_argument('--input', action='store', required=True,
                          help='binary excel file', dest='input')
    parser_e.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_e.add_argument('--output', action='store', required=True,
                          help='directory containing results', dest='output')

    parser_t = subparsers.add_parser('transpose', help='transpose table')
    parser_t.add_argument('--input', action='store', required=True,
                          help='text table file', dest='input')
    parser_t.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_t.add_argument('--output', action='store', required=True,
                          help='transposed table', dest='output')

    parser_f = subparsers.add_parser('fillna', help='fill NA values')
    parser_f.add_argument('--input', action='store', required=True,
                          help='text table file', metavar='<FILE>', dest='input')
    parser_f.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_f.add_argument('--nav', action='store', choices=['NA', 'NaN', 'Null', '0', 'blank'], required=True,
                          help='strings to recognize as NA', dest='nav')
    parser_f.add_argument('--rmp', action='store', type=float,
                          help='drop row with more than specified percentage of NA values',
                          metavar='<INT>', dest='rmp')
    parser_f.add_argument('--method', action='store',
                          choices=['min', 'otMin', 'halfMin', 'colMin', 'colMean', 'missForest', 'kNN', 'PPCA', 'BPCA',
                                   'nipals'], required=True,
                          help='method for filling NA values', dest='method')
    parser_f.add_argument('--r_interpreter', action='store', required=True,
                          help='path of R interpreter', metavar='<FILE>', dest='r_interpreter')
    parser_f.add_argument('--r_script', action='store', required=True,
                          help='requried R script', metavar='<FILE>', dest='r_script')
    parser_f.add_argument('--output', action='store', required=True,
                          help='NA filled table', metavar='<FILE>', dest='output')

    parser_s = subparsers.add_parser('standard', help='standardize matrix')
    parser_s.add_argument('--input', action='store', required=True,
                          help='text table file', metavar='<FILE>', dest='input')
    parser_s.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_s.add_argument('--observe', action='store', choices=['none', 'sum', 'median', 'mean', 'refcol', 'refrow'], required=True,
                          help='method for strandardizing observations', dest='observe')
    parser_s.add_argument('--guard', action='store',
                          help='specified index name when selecting refcol or refrow', metavar='<STR>', dest='guard')
    parser_s.add_argument('--feature', action='store', choices=['none', 'standard_scale', 'minmax_scale', 'maxabs_scale', 'normalize'], required=True,
                          help='sklearn method for strandardizing matrix', dest='feature')
    parser_s.add_argument('--output', action='store', required=True,
                          help='transposed table', metavar='<FILE>', dest='output')

    parser_c = subparsers.add_parser('filter', help='filter tables')
    parser_c.add_argument('--input', action='store', required=True,
                          help='text table file', metavar='<FILE>', dest='input')
    parser_c.add_argument('--sep', action='store', choices=['tab', 'comma', 'semicolon', 'space'], required=True,
                          help='separator type', dest='sep')
    parser_c.add_argument('--table_head', action='store', required=True, type=bool,
                          help='if contain a header', dest='table_head')
    parser_c.add_argument('--contain', action='store', required=True, type=bool,
                          help='if contain', dest='contain')
    parser_c.add_argument('--lower', action='store', required=True, type=bool,
                          help='if lower', dest='lower')
    parser_c.add_argument('--number', action='store', required=True,
                          help='if contain a header', dest='number')
    parser_c.add_argument('--info', action='store', required=True,
                          help='filter info', dest='info')
    parser_c.add_argument('--output', action='store', required=True,
                          help='filtered table', dest='output')


    logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

    if len(sys.argv) > 1 and sys.argv[1] not in ['-h', '--help']:
        if sys.argv[1] in ['combine', 'excel', 'transpose', 'fillna', 'standard', 'filter']:
            args = parser.parse_args()
        else:
            logging.error('unrecognized command ({})'.format(sys.argv[1]))
            sys.exit(-2)
    else:
        parser.print_help()
        sys.exit(-1)

    inst = TableKit(args)
    inst.run()
