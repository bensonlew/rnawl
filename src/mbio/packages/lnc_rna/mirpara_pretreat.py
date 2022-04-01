# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import re
import os

parser = OptionParser(description='Conduct id conversion for input FASTA file')
parser.add_option('-i', '--input', dest='input', help='input raw lncRNA FASTA file')
parser.add_option('-m', '--map', dest='map', help='output relation map file about new id and old id')
parser.add_option('-o', '--output', dest='output', help='output pretreated lncRNA FASTA file')
(opts, args) = parser.parse_args()

def main(file_in, map_out, file_out):
    print 'INFO: start processing {}'.format(file_in)
    map_lines = list()
    out_lines = list()
    for line in open(file_in):
        if line[0] == '>':
            old_id = line.split()[0][1:]
            new_id = old_id.replace('_', '')
            map_lines.append('{}\t{}\n'.format(new_id, old_id))
            out_lines.append(re.sub(old_id, new_id, line, count=1))
        else:
            out_lines.append(line)
    open(map_out, 'w').writelines(map_lines)
    open(file_out, 'w').writelines(out_lines)
    if os.path.getsize(map_out) > 0:
        print 'INFO: succeed in exporting {}'.format(map_out)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.map and opts.output:
        main(opts.input, opts.map, opts.output)
    else:
        parser.print_help()