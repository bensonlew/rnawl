# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Generate sequence id list file from specified input')
parser.add_option('--input', dest='input', type=str, help='input single annotation result file')
parser.add_option('--database', dest='database', type=str, help='input database of annotation')
parser.add_option('--output', dest='output', type=str, help='output id list file')
(opts, args) = parser.parse_args()

def main(file_in, database, file_out):
    print 'INFO: start reading {}'.format(file_in)
    open(file_out, 'w').writelines('{}\n'.format(i) for i in process(file_in, database))
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

def process(file_in, database):
    if database == 'nr':
        ids = set(pd.read_table(file_in)['Query-Name'])
    elif database == 'swissprot':
        ids =  set(pd.read_table(file_in)['Query-Name'])
    elif database == 'cog':
        ids = set()
        ids.update(*pd.read_table(file_in)['Seqs'].apply(lambda s: set(s.split(';'))))
    elif database == 'kegg':
        ids = set(line.strip().split('\t')[0] for n, line in enumerate(open(file_in)) if n)
    elif database == 'pfam':
        ids = set(pd.read_table(file_in)['Seq_id'])
    elif database == 'go':
        ids = set(line.strip().split('\t')[0] for line in open(file_in))
    return ids

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['input', 'database', 'output'])):
        main(opts.input, opts.database, opts.output)
    else:
        parser.print_help()
