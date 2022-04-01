# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import os
import glob
from Bio import motifs
from Bio import SeqIO
# from Bio import Alphabet
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
import re
from biocluster.core.exceptions import OptionError


class Seqlogo(object):
    def __init__(self, infile, m_type, unit, s_type):
        self.infile = infile
        self.unit = unit
        self.m_type = m_type.lower()
        self.s_type = s_type

    def parse_seq(self):
        seq_list = list()
        for each in SeqIO.parse(self.infile, 'fasta'):
            sequence = str(each.seq).upper()
            seq_list.append(Seq(sequence))
        motif = motifs.create(seq_list)
        self.draw_weblogo(seq_type=self.s_type, seq_id='motif', motif=motif, unit=self.unit)
        return motif

    def parse_matrix(self):
        try:
            with open(self.infile) as handle:
                record = motifs.parse(handle, self.m_type, strict=False)
        except:
            if self.m_type == 'meme':
                self.m_type = 'minimal'
                with open(self.infile) as handle:
                    record = motifs.parse(handle, 'minimal', strict=False)
            elif self.m_type == 'pfm':
                self.m_type = 'jaspar'
                with open(self.infile) as handle:
                    record = motifs.parse(handle, 'jaspar', strict=False)
            else:
                raise ValueError('Unexpected format of matrix.')
        i = 1
        for each in record:
            if self.m_type == 'transfac':
                seq_id = each['AC']
            elif self.m_type in ['meme', 'minimal', 'pfm-four-columns', 'xms']:
                seq_id = each.name
            elif self.m_type == 'jaspar':
                seq_id = each.matrix_id
            else:
                seq_id = 'motif_' + str(i)
            self.draw_weblogo(seq_type=self.s_type, seq_id=seq_id, motif=each, unit=self.unit)
            i += 1

    # def process_transfac(self, transfac):
    #     result = transfac.read()
    #     transfac_path = 'jaspar_transfac.dat'
    #     with open(transfac_path, "w") as out:
    #         for i in result.content.split('\n'):
    #             item = re.split('\t| ', i, maxsplit=1)
    #             if len(item) > 1:
    #                 out.write(item[0] + '\t\t' + item[1] + '\n')
    #             else:
    #                 out.write(item[0] + '\n')
    #     with open(transfac) as handle:
    #         record = motifs.parse(handle, 'transfac')
    #     # os.remove(transfac)
    #     return record

    # def generate_consensus(self, m_type, motif):
    #     if m_type.lower() == 'pfm':
    #         consensus = motif.counts.consensus
    #     else:
    #         consensus = motif.consensus
    #     seq_type = self.verify_type(consensus=consensus)
    #     return seq_type
    #
    # def verify_type(self, consensus):
    #     type_list = [IUPAC.unambiguous_dna, IUPAC.unambiguous_rna, IUPAC.protein, IUPAC.ambiguous_dna, IUPAC.ambiguous_rna]
    #     for i in range(5):
    #         seq = Seq(str(consensus).upper(), type_list[i])
    #         test = Alphabet._verify_alphabet(seq)
    #         if str(test) == 'True' and i in [0, 3]:
    #             seq_type = 'dna'
    #             break
    #         elif str(test) == 'True' and i in [1, 4]:
    #             seq_type = 'rna'
    #             break
    #         elif str(test) == 'True' and i == 2:
    #             seq_type = 'protein'
    #             break
    #         elif str(test) == 'False' and i == 4:
    #             print('Unknown motif type of {}'.format(consensus))
    #             seq_type = ''
    #             break
    #         else:
    #             continue
    #     return seq_type

    def draw_weblogo(self, seq_type, seq_id, motif, unit='both'):
        if seq_type == 'nt':
            color = 'color_classic'
        elif seq_type == 'aa':
            color = 'color_chemistry'
        else:
            color = 'color_auto'
        if unit == 'both':
            motif.weblogo('{}_bits.pdf'.format(seq_id), format='pdf', show_fineprint=False, show_errorbars=False,
                          unit_name='bits', stack_width='large', alphabet='alphabet_auto', color_scheme=color,
                          scale_width=False)
            unit = 'probability'
        motif.weblogo('{}_{}.pdf'.format(seq_id, unit), format='pdf', show_fineprint=False, show_errorbars=False,
                      unit_name=unit, stack_width='large', alphabet='alphabet_auto', color_scheme=color, scale_width=False)

    def run_seqlogo(self):
        if self.m_type in ['none', 'sites']:
            self.parse_seq()
        else:
            self.parse_matrix()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='A script to draw sequence logo for motifs')
    parser.add_argument('-i', type=str, help='A matrix file or a sequence file')
    parser.add_argument('-m', type=str, help='Type of matrix. If a sequence file, the type of matrix is none. Otherwise, choose one from PFM, MEME, JASPAR, TRANSFAC.')
    parser.add_argument('-u', type=str, help='Units of y-axis')
    parser.add_argument('-s', type=str, help='Sequence type, one of nt(RNA/DNA), ')
    args = parser.parse_args()
    conversion = Seqlogo(args.i, args.m, args.u, args.s)
    conversion.run_seqlogo()
