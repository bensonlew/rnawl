# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os

parser = OptionParser(description='Change chrom from input BED by input chromosome tabular MAP')
parser.add_option('-i', '--input', dest='input', help='input raw BED file')
parser.add_option('-m', '--map', dest='map', help='input relation map file about old chrom and new chrom')
parser.add_option('-o', '--output', dest='output', help='output processed BED file')
(opts, args) = parser.parse_args()

def main(file_in, map, file_out):
    print 'INFO: start reading {}'.format(map)
    chrom_map_dict = dict([line.strip().split('\t') for line in open(map)])
    print 'INFO: start reading {}'.format(file_in)
    lines = list()
    for line in open(file_in):
        items = line.split('\t')
        if len(items) > 0 and items[0] in chrom_map_dict:
            lines.append('\t'.join([chrom_map_dict[items[0]]] + items[1:]))
    else:
        open(file_out, 'w').writelines(lines)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.map and opts.output:
        main(opts.input, opts.map, opts.output)
    else:
        parser.print_help()
