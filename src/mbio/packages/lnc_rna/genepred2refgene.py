# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Convert genePred file to refGene file')
parser.add_option('-i', '--input', dest='input', help='Input genePred file')
parser.add_option('-o', '--output', dest='output', help='Output refGene file')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start reading {}'.format(file_in)
    with open(file_out, 'w') as f:
        for line in open(file_in):
            items = line.strip().split('\t')
            if len(items) == 15 and items[2] in ['+', '-']:
                f.write('1\t{}'.format(line))
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()