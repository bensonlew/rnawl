# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import re
import os

parser = OptionParser(description='Restore transcript ids of miRPara result file')
parser.add_option('-i', '--input', dest='input', help='input raw miRPara out file')
parser.add_option('-m', '--map', dest='map', help='input relation map file about new id and old id')
parser.add_option('-o', '--output', dest='output', help='output processed miRPara out file')
(opts, args) = parser.parse_args()

def main(file_in, map_in, file_out):
    print 'INFO: start processing {}'.format(file_in)
    map_dict = dict([line.strip().split('\t') for line in open(map_in)])
    lines = list()
    for line in open(file_in):
        if line[0] != '#':
            items = line.split('\t')
            new_id = items[0].split(':')[0]
            old_id = map_dict[new_id]
            lines.append(line.replace(new_id, old_id))
    open(file_out, 'w').writelines(lines)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.map and opts.output:
        main(opts.input, opts.map, opts.output)
    else:
        parser.print_help()