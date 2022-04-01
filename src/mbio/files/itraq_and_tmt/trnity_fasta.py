# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.iofile import File
from collections import defaultdict
import re
import subprocess
from biocluster.config import Config
import os
from biocluster.core.exceptions import FileError
from Bio import SeqIO
from mbio.files.denovo_rna_v2.fasta import FastaFile
import unittest

class TrinityFasta(FastaFile):
    """
    定义Trinity Fasta文件
    """

    def __init__(self):
        super(TrinityFasta, self).__init__()
        self._is_unigene = False
        self.trans2gene = {}
        self.unigene_trans = {}
        self.seq = {}

    def set_gene2tran(self, gene2trans):
        """
        设置基因转录本对应关系
        """
        with open(gene2trans, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                if line[1]:
                    self.trans2gene[line[1]] = line[0]
                else:
                    self.set_error("gene2trans至少包含两列", code="42501701")
        return self.gene2trans

    def get_unigene(self, output, marker_prefix):
        """
        提取unigene序列文件,和对应关系标记，长度等信息
        """
        if self._is_unigene == True:
            self.set_error("文件必须是转录本序列", code="42501702")
        elif not self.trans2gene:
            self.set_error("需使用gene2trans设置基因转录本对应关系文件", code="42501703")
        else:
            trans_len = self.get_contig_len()
            for trans in trans_len.keys():
                gene = self.trans2gene[trans]
                if self.unigene_trans.has_key(gene):
                    if self.unigene_trans[gene]['len'] < trans_len[trans]:
                        self.unigene_trans[gene]['trans'] = trans
                        self.unigene_trans[gene]['len'] = trans_len[trans]
                    else:
                        pass
                else:
                    self.unigene_trans[gene]['trans'] = trans
                    self.unigene_trans[gene]['len'] = trans

            seq_records = SeqIO.parse(self.path, 'fasta')
            for seq_record in seq_records:
                seq_seq = seq_record.seq
                seq_name = seq_record.name
                self.seq[seq_name] = seq_seq

            with open(output, 'w') as unigene:
                for gene in self.unigene_trans.keys():
                    trans = self.unigene_trans[gene]['trans']
                    unigene.write('>{}\n{}\n'.format(gene, self.seq[trans]))


            with open(marker_prefix + '_t2g', 'w') as mar1, \
                 open(marker_prefix + '_t2g2u', 'w') as mar2, \
                 open(marker_prefix + '_t2g2u2l', 'w') as mar3:

                for trans, gene in self.trans2gene.items():
                    marker = 'no'
                    if self.unigene_trans[gene]['trans'] == trans:
                        marker = 'yes'
                    mar1.write('{}\t{}\n'.format(trans, gene))
                    mar2.write('{}\t{}\t{}\n'.format(trans, gene, marker))
                    mar3.write('{}\t{}\t{}\t{}\n'.format(trans, gene, marker, trans_len[trans]))



class TestFunction(unittest.TestCase):
    """
    This is test for the trinityfasta
    """
    def test(self):
        fasta = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/trinity_out_dir/Trinity.fasta'
        g2t = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/trinity_out_dir/Trinity.fasta.gene_trans_map'

        fa = TrinityFasta()
        fa.prop['pathway'] = fasta
        fa.set_gene2tran(g2t)
        fa.get_unigene(fasta + '.unigene.fa', fasta)

if __name__ == '__main__':
    unittest.main()
