# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
from Bio import SeqIO
import os

parser = OptionParser(description='Split genome FASTA to chromosome BED')
parser.add_option('-i', '--input', dest='input', help='Input genome FASTA file')
parser.add_option('-o', '--output', dest='output', help='Output directory for chromosome BED files')
(opts, args) = parser.parse_args()

def main(file_in, dir_out):
    print 'INFO: start processing {}'.format(file_in)
    for n, seq_record in enumerate(SeqIO.parse(file_in, 'fasta')):
        bed_out = os.path.join(dir_out, 'chr{}.bed'.format(seq_record.id))
        open(bed_out, 'w').write('{}\t0\t{}\n'.format(seq_record.id, len(seq_record)))
        if os.path.getsize(bed_out) > 0:
            print 'INFO: succeed in exporting {}'.format(bed_out)

if __name__ == '__main__':
    if opts.input and opts.output:
        main(opts.input, opts.output)
    else:
        parser.print_help()