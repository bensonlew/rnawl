# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
from Bio import SeqIO
import os

parser = OptionParser(description='Filter FASTA file by input id list')
parser.add_option('-i', '--input', dest='input', help='input raw FASTA file')
parser.add_option('-l', '--list', dest='list', help='input file containing id that needs retaining')
parser.add_option('-o', '--output', dest='output', help='output filtered FASTA file')
(opts, args) = parser.parse_args()

def main(file_in, id_list, file_out):
    print 'INFO: start reading {}'.format(id_list)
    ids = {line.strip() for line in open(id_list)}
    print 'INFO: start reading {}'.format(file_in)
    SeqIO.write([seq_record for seq_record in SeqIO.parse(file_in, 'fasta') if seq_record.id in ids], file_out, 'fasta')
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.list and opts.output:
        main(opts.input, opts.list, opts.output)
    else:
        parser.print_help()
