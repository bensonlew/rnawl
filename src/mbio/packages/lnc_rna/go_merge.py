# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os

parser = OptionParser(description='Convert blast2go annotation table to go list table')
parser.add_option('-i', '--input', dest='input', help='Input blast2go annotation file')
parser.add_option('-o', '--output', dest='output', help='Output go list file')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start processing {}'.format(file_in)
    dct = dict()
    for n, line in enumerate(open(file_in)):
        items = line.split('\t')
        if len(items) < 3:
            print 'WARNING: drop {} at line {}'.format(line, n + 1)
            continue
        seq_name = items[0]
        go_list = items[1].split(';')
        if dct.has_key(seq_name):
            dct[seq_name].extend(go_list)
        else:
            dct[seq_name] = go_list
    with open(file_out, 'w') as f:
        for k, v in dct.iteritems():
            f.write('{}\t{}\n'.format(k, ';'.join(set(v))))
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()