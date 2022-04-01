# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Generate gene PFAM DOMAIN file')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--domain', dest='domain', help='input transcript PFAM DOMAIN file')
parser.add_option('--output', dest='output', help='output gene PFAM DOMAIN file')
(opts, args) = parser.parse_args()

def main(t2g_tsv, txpt_tsv, gene_tsv):
    print 'INFO: start reading {}'.format(t2g_tsv)
    t2g_dict = dict(line.strip().split('\t')[:2] for line in open(t2g_tsv))
    print 'INFO: start reading {}'.format(txpt_tsv)
    df = pd.read_table(txpt_tsv, header=0).query('Seq_id in @t2g_dict')
    df['Seq_id'] = df['Seq_id'].apply(lambda t: t2g_dict[t])
    df.to_csv(gene_tsv, sep='\t', index=None)
    if os.path.getsize(gene_tsv) > 0:
        print 'INFO: succeed in exporting {}'.format(gene_tsv)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['t2g', 'domain', 'output'])):
        main(opts.t2g, opts.domain, opts.output)
    else:
        parser.print_help()
