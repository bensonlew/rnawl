# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os

parser = OptionParser(description='Filter BED file by input name list')
parser.add_option('-i', '--input', dest='input', help='input raw BED file')
parser.add_option('-n', '--name', dest='name', help='input file containing name that needs retaining')
parser.add_option('-o', '--output', dest='output', help='output filtered BED file')
(opts, args) = parser.parse_args()

def main(file_in, name, file_out):
    print 'INFO: start reading {}'.format(name)
    names = [line.strip() for line in open(name)]
    print 'INFO: start reading {}'.format(file_in)
    lines = list()
    for line in open(file_in):
        items = line.split('\t')
        if len(items) > 3 and items[3] in names:
            lines.append(line)
    else:
        open(file_out, 'w').writelines(lines)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.name and opts.output:
        main(opts.input, opts.name, opts.output)
    else:
        parser.print_help()
