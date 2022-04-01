# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng, liubinxu'

from optparse import OptionParser
import os

parser = OptionParser(description='Generate gene reactome_ids annotation file')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--disgenet_ids', dest='disgenet_ids', help='input transcript disgenet_ids annotation file')
parser.add_option('--output', dest='output', help='output gene disgenet_ids annotation file')
(opts, args) = parser.parse_args()

def main(t2g_tsv, txpt_b2g, gene_b2g):
    print 'INFO: start reading {}'.format(t2g_tsv)
    t2g_dict = dict(line.strip().split('\t') for line in open(t2g_tsv))
    print 'INFO: start reading {}'.format(txpt_b2g)
    lines = list()
    for n, line in enumerate(open(txpt_b2g)):
        items = line.strip().split('\t')
        if len(items) < 2:
            print 'WARN: drop line {}\n{}'.format(n + 1, line)
            continue
        if items[0] in t2g_dict:
            lines.append('{}\n'.format('\t'.join([t2g_dict[items[0]]] + items[1:])))
    else:
        open(gene_b2g, 'w').writelines(lines)
        if os.path.getsize(gene_b2g) > 0:
            print 'INFO: succeed in exporting {}'.format(gene_b2g)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['t2g', 'disgenet_ids', 'output'])):
        main(opts.t2g, opts.disgenet_ids, opts.output)
    else:
        parser.print_help()
