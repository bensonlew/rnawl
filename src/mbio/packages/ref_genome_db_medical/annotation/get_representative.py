# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Generate relationship file between longest transcript and gene')
parser.add_option('--input', dest='input', type=str, help='input T2G2R2L2P file')
parser.add_option('--output', dest='output', type=str, help='output longest T2G file')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start reading {}'.format(file_in)
    df = pd.read_table(file_in, header=None)
    df[df[2] == 'yes'].reindex([0, 1], axis=1).to_csv(file_out, sep='\t', header=False, index=False)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 2, ['input', 'output'])):
        main(opts.input, opts.output)
    else:
        parser.print_help()
