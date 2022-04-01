# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
from Bio import SeqIO
import os

parser = OptionParser(description='Split FASTA for CMSCAN')
parser.add_option('--input', dest='input', help='input FASTA file')
parser.add_option('--number', dest='number', type=int, help='number of the split FASTA file')
parser.add_option('--output', dest='output', help='output directory containing results')
(opts, args) = parser.parse_args()

def main(file_in, number, dir_out):
    print 'INFO: start reading {}'.format(file_in)
    seq_records = list(SeqIO.parse(file_in, 'fasta'))
    if len(seq_records) < number:
        SeqIO.write(seq_records, os.path.join(dir_out, 'lncrna.fasta'), 'fasta')
    else:
        step = len(seq_records) / number
        for n, i in enumerate(range(0, len(seq_records), step)):
            if n + 1 < number:
                SeqIO.write(seq_records[i: i + step], os.path.join(dir_out, 'lncrna.{}.fasta'.format(n)), 'fasta')
            elif n + 1 == number:
                SeqIO.write(seq_records[i:], os.path.join(dir_out, 'lncrna.{}.fasta'.format(n)), 'fasta')
    for basename in os.listdir(dir_out):
        if os.path.getsize(os.path.join(dir_out, basename)) > 0:
            print 'INFO: succeed in exporting {}'.format(os.path.join(dir_out, basename))

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 3, ['input', 'number', 'output'])):
        main(opts.input, opts.number, opts.output)
    else:
        parser.print_help()
