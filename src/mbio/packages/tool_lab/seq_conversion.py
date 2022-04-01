# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from biocluster.core.exceptions import OptionError


class Conversion(object):
    def __init__(self, infile, dna_type, output):
        self.infile = infile
        self.dna_dict = dict()
        self.dna_product = dict()
        self.rna_dict = dict()
        self.rna_product = dict()
        self.seq_type = [IUPAC.unambiguous_dna, IUPAC.unambiguous_rna, IUPAC.ambiguous_dna, IUPAC.ambiguous_rna]
        self.dna_type = dna_type
        self.output = output

    def get_order(self):
        order_list = list()
        with open(self.infile, 'r') as fa:
            for line in fa:
                if line.startswith('>'):
                    order_list.append(line.strip())
        return order_list

    def process_file(self):
        record = SeqIO.parse(self.infile, 'fasta')
        for each in record:
            self.process_each(sequence=each.seq,
                              title=each.id)

    def process_each(self, sequence, title):
        for i in range(4):
            seq = Seq(str(sequence).upper(), self.seq_type[i])
            test = Alphabet._verify_alphabet(seq)
            if str(test) == 'True' and i in [0, 2]:
                self.dna_dict[title] = sequence
                break
            elif str(test) == 'True' and i in [1, 3]:
                self.rna_dict[title] = sequence
                break
            elif str(test) == 'False' and i == 3:
                print('Unknown sequence type of {}:{}'.format(title, str(sequence)))
                break
            else:
                continue

    def dna2rna(self):
        if self.dna_type == 'coding':
            for k, v in self.dna_dict.items():
                self.dna_product[k] = v.transcribe()
        if self.dna_type == 'template':
            for k, v in self.dna_dict.items():
                self.dna_product[k] = v.reverse_complement().transcribe()

    def rna2dna(self):
        for k, v in self.rna_dict.items():
            self.rna_product[k] = v.back_transcribe()

    def merge(self):
        product = dict(self.dna_product.items() + self.rna_product.items())
        order_list = self.get_order()
        with open(self.output, 'w') as o:
            for each in order_list:
                o.write(each + '\n')
                key = each.split(' ')[0]
                v = product.get(key[1:])
                if v:
                    value = str(v)
                    while len(value) > 60:
                        o.write(value[0:60] + '\n')
                        value = value[60:len(value)]
                    o.write(value + '\n')
                else:
                    o.write('Improper sequence type\n')

    def run_conversion(self):
        self.process_file()
        self.dna2rna()
        self.rna2dna()
        self.merge()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='A script to mimic DNA/RNA (reverse-)transcription.')
    parser.add_argument('-i', type=str, help='Sequence file in fasta format')
    parser.add_argument('-t', type=str, help='Type of DNA strand. Choose one from coding/template.')
    parser.add_argument('-o', type=str, help='output path')
    args = parser.parse_args()
    conversion = Conversion(args.i, args.t, args.o)
    conversion.run_conversion()




