# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Process results of lncRNA location conservation analysis')
parser.add_option('-i', '--input', dest='input', help='input BWTOOL summary out file')
parser.add_option('-r', '--retain', dest='retain', help='input file containing name that needs retaining')
parser.add_option('-t', '--type', dest='type', help='input file containing name and its corresponding type')
parser.add_option('-o', '--output', dest='output', help='output tabular file')
(opts, args) = parser.parse_args()

def main(file_in, retain_in, type_in, file_out):
    print 'INFO: start reading {}'.format(retain_in)
    names = [line.strip() for line in open(retain_in)]
    print 'INFO: start reading {}'.format(file_in)
    df = pd.read_table(file_in)
    df = df.query('name in @names')
    df = df.reindex(
        ['name', '#chrom', 'start', 'end', 'size', 'num_data', 'min', 'max', 'mean', 'median', 'sum'], axis=1
    )
    df = df.rename(columns={'name': 'lncrna_id', '#chrom': 'chromosome', 'num_data': 'number'})
    df = df.merge(pd.read_table(type_in, names=['lncrna_id', 'lncrna_type'], usecols=[0, 3]), on='lncrna_id', how='left')
    df.to_csv(file_out, sep='\t', index=False)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.retain and opts.type and opts.output:
        main(opts.input, opts.retain, opts.type, opts.output)
    else:
        parser.print_help()
