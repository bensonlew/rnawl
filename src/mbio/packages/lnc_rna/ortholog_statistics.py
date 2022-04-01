# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Conduct statistics of lncRNA ortholog analysis results')
parser.add_option('-f', '--forward', dest='forward', help='input forward blastn outfmt6 file')
parser.add_option('-r', '--reverse', dest='reverse', help='input reverse blastn outfmt6 file')
parser.add_option('-i', '--identity', dest='identity', type=float, help='the lower limit of identity')
parser.add_option('-t', '--type', dest='type', help='input file containing name and its corresponding type')
parser.add_option('-o', '--output', dest='output', help='output tabular file')
(opts, args) = parser.parse_args()

def main(forward_tblout, reverse_tblout, identity, type_in, file_out):
    print 'INFO: start processing {}'.format(reverse_tblout)
    reverse_qs_pairs = [line.strip().split('\t')[:2][::-1] for line in open(reverse_tblout)]
    print 'INFO: start processing {}'.format(forward_tblout)
    df = pd.read_table(forward_tblout, names=[
        'query_id', 'subject_id', 'identity', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end',
        's_start', 's_end', 'e_value', 'bit_score'
    ])
    df = df.reindex([i for i in range(df.shape[0]) if list(df.iloc[i][:2]) in reverse_qs_pairs])
    df = df[df['identity'] >= identity]
    df = df.reindex(['query_id', 'subject_id', 'e_value', 'identity'], axis=1)
    df = df.rename(columns={'query_id': 'lncrna_id', 'subject_id': 'orthologous'})
    df = df.merge(pd.read_table(type_in, names=['lncrna_id', 'lncrna_type'], usecols=[0, 3]), on='lncrna_id', how='left')
    df.to_csv(file_out, sep='\t', index=False)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 5, ['forward', 'reverse', 'identity', 'type', 'output'])):
        main(opts.forward, opts.reverse, opts.identity, opts.type, opts.output)
    else:
        parser.print_help()
