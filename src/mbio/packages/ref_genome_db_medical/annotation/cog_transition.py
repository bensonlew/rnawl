# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng, liubinxu'

from optparse import OptionParser
import os

parser = OptionParser(description='Generate gene cog annotation file')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--cog', dest='cog', help='input transcript cog annotation file')
parser.add_option('--output', dest='output', help='output gene cog annotation file')
(opts, args) = parser.parse_args()

def main(t2g_tsv, txpt_cog, gene_cog):
    print 'INFO: start reading {}'.format(t2g_tsv)
    t2g_dict = dict(line.strip().split('\t') for line in open(t2g_tsv))
    print 'INFO: start reading {}'.format(txpt_cog)
    lines = list()
    for n, line in enumerate(open(txpt_cog)):
        items = line.strip().split('\t')
        if len(items) < 2:
            print 'WARN: drop line {}\n{}'.format(n + 1, line)
            continue
        if items[0] in t2g_dict:
            lines.append('{}\n'.format('\t'.join([t2g_dict[items[0]]] + items[1:])))
        elif line[0].startswith("#"):
            lines.append(line)
    else:
        open(gene_cog, 'w').writelines(lines)
        if os.path.getsize(gene_cog) > 0:
            print 'INFO: succeed in exporting {}'.format(gene_cog)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['t2g', 'cog', 'output'])):
        main(opts.t2g, opts.cog, opts.output)
    else:
        parser.print_help()
