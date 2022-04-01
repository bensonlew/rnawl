# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Filter rMATS output event file by provided gene id')
parser.add_option('-i', '--input', dest='input', help='input raw event file')
parser.add_option('-g', '--gene', dest='gene', help='input gene id that needs retaining')
parser.add_option('-o', '--output', dest='output', help='output filtered event file')
(opts, args) = parser.parse_args()

def main(file_in, gene, file_out):
    print 'INFO: start reading {}'.format(file_in)
    with open(file_out, 'w') as f:
        for n, line in enumerate(open(file_in)):
            if n == 0:
                f.write(line)
            elif len(line.split('\t')) > 1 and line.split('\t')[1] == gene:
                f.write(line)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.gene and opts.output:
        main(opts.input, opts.gene, opts.output)
    else:
        parser.print_help()
