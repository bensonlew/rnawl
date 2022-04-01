# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Get tabular file containing correspondence bewteen id, name and description')
parser.add_option('-i', '--input', dest='input', help='Input BIOMART file')
parser.add_option('-t', '--type', dest='type', help='The type of BIOMART file')
parser.add_option('-o', '--output', dest='output', help='Output tabular file')
(opts, args) = parser.parse_args()

def main(file_in, des_type, file_out):
    print 'INFO: start reading {}'.format(file_in)
    df = pd.DataFrame([], columns=['gene_id', 'symbol', 'des'])
    lines = list()
    for n, line in enumerate(open(file_in)):
        if line.strip() != '':
            lines.append(line_to_list(line, des_type))
    df = pd.DataFrame(lines, columns=['gene_id', 'symbol', 'des'])
    df = df.drop_duplicates()
    df.to_csv(file_out, index=False, sep='\t', header=True)
    if os.path.getsize(file_out):
        print 'INFO: succeed in exporting {}'.format(file_out)

def line_to_list(line, des_type=None):
    items = line.strip().split('\t')
    gene_id = items[0] if items[0] else '-'
    if not des_type:
        return
    elif des_type == 'type1':
        symbol = items[2] if items[2] else '-'
        des = items[7] if items[7] else '-'
    elif des_type == 'type2':
        symbol = items[2] if items[2] else '-'
        des = items[5] if items[5] else '-'
    elif des_type == 'type3':
        symbol = '-'
        des = items[3] if items[3] else '-'
    return [gene_id, symbol, des]

if __name__ == '__main__':
    if opts.input and opts.type and opts.output:
        main(opts.input, opts.type, opts.output)
    else:
        parser.print_help()