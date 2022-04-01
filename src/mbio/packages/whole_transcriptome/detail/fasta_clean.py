# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import re

from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    p = re.compile(r'\(\S\)')
    sequences = list()
    for record in SeqIO.parse(args.ifasta, 'fasta'):
        record.id = re.sub(p, '', record.id)
        record.name = re.sub(p, '', record.name)
        record.description = re.sub(p, '', record.description)
        sequences.append(record)
    SeqIO.write(sequences, args.ofasta, 'fasta')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Clean bedtools output FASTA file')
    parser.add_argument('-i', action='store', required=True,
                        help='original FASTA', metavar='<FILE>', dest='ifasta')
    parser.add_argument('-o', action='store', required=True,
                        help='formated FASTA', metavar='<FILE>', dest='ofasta')

    args = parser.parse_args()

    main(args)
