# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os

parser = OptionParser(description='Export relationship file between id and GO terms from blast2go result')
parser.add_option('-i', '--input', dest='input', help='input blast2go annotation file')
parser.add_option('-o', '--output', dest='output', help='output relationship file between id and GO terms')
(opts, args) = parser.parse_args()

def main(file_in, file_out):
    print 'INFO: start processing {}'.format(file_in)
    dct = dict()
    for n, line in enumerate(open(file_in)):
        items = line.strip().split('\t')
        if len(items) < 2:
            print 'WARN: drop line {}\n{}'.format(n + 1, line)
            continue
        seq_id = items[0]
        go_list = items[1].split(';')
        if dct.has_key(seq_id):
            dct[seq_id].extend(go_list)
        else:
            dct[seq_id] = go_list
    else:
        open(file_out, 'w').writelines('{}\t{}\n'.format(k, ';'.join(set(v))) for k, v in dct.items())
        if os.path.getsize(file_out) > 0:
            print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()
