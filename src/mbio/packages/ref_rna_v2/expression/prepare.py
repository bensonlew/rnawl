# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import numpy as np
from Bio import SeqIO
import os

parser = OptionParser(description='Prepare transcripts FASTA file for quantification')
parser.add_option('-i', dest='input', type=str, help='input raw transcripts FASTA file')
parser.add_option('-m', dest='txpt2gene', type=str, help='input raw file containing the correspondent relationship between transcript and gene')
parser.add_option('-t', dest='t2g', type=str, help='output clean T2G relationship file')
parser.add_option('-g', dest='g2t', type=str, help='output clean G2T relationship file')
parser.add_option('-o', dest='output', type=str, help='output clean transcripts FASTA file')
(opts, args) = parser.parse_args()

def main(fasta_in, map, t2g, g2t, fasta_out):
    print 'INFO: start reading {}'.format(map)
    map_dict = dict([line.strip().split('\t') for line in open(map) if len(line.strip().split('\t')) == 2])
    print 'INFO: start reading {}'.format(fasta_in)
    record_dict = dict([(seq_record.id, seq_record) for seq_record in SeqIO.parse(fasta_in, 'fasta')])
    txpt_arr = np.intersect1d(map_dict.keys(), record_dict.keys())
    seq_records = list()
    with open(t2g, 'w') as t2g_handle, open(g2t, 'w') as g2t_handle:
        for txpt in txpt_arr:
            t2g_handle.write('{}\t{}\n'.format(txpt, map_dict[txpt]))
            g2t_handle.write('{}\t{}\n'.format(map_dict[txpt], txpt))
            seq_records.append(record_dict[txpt])
    SeqIO.write(seq_records, fasta_out, 'fasta')
    if os.path.getsize(t2g) > 0:
        print 'INFO: succeed in exporting {}'.format(t2g)
    if os.path.getsize(g2t) > 0:
        print 'INFO: succeed in exporting {}'.format(g2t)
    if os.path.getsize(fasta_out) > 0:
        print 'INFO: succeed in exporting {}'.format(fasta_out)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 5, ['input', 'txpt2gene', 'g2t', 't2g', 'output'])):
        main(opts.input, opts.txpt2gene, opts.t2g, opts.g2t, opts.output)
    else:
        parser.print_help()
