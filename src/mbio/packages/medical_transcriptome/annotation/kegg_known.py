# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import os

parser = OptionParser(description='Generate KNOWN KO file for KEGG gene annotation')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--ko', dest='ko', type=str, help='input raw KNOWN KO file')
parser.add_option('--output', dest='output', type=str, help='output KNOWN KO file for gene annotation')
(opts, args) = parser.parse_args()

def main(t2g_tsv, ko_in, ko_out):
    print 'INFO: start reading {}'.format(t2g_tsv)
    t2g_dict = dict(line.strip().split('\t') for line in open(t2g_tsv))
    print 'INFO: start reading {}'.format(ko_in)
    lines = list()
    for line in open(ko_in):
        items = line.strip().split('\t')
        if len(items) > 2 and items[1] in t2g_dict:
            items[1] = items[0]
            lines.append('{}\n'.format('\t'.join(items)))
    else:
        open(ko_out, 'w').writelines(lines)
        if os.path.getsize(ko_out) > 0:
            print 'INFO: succeed in exporting {}'.format(ko_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['t2g', 'ko', 'output'])):
        main(opts.t2g, opts.ko, opts.output)
    else:
        parser.print_help()
