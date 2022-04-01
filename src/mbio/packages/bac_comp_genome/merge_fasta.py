# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os,re
import pandas as pd
import argparse
from Bio import SeqIO


class ortholog(object):
    """
    根据泛基因组的cluster表，提取单拷贝的基因蛋白或核苷酸序列
    """
    def get_merge_fasta(self, dir, sample_list, out):
        list1 = sample_list.split(",")
        files = os.listdir(dir)
        with open (out, "w") as g:
            for sample in list1:
                seqs = ''
                for file in files:
                    seq = self.get_gene_id(dir + "/" + file)
                    if sample in seq:
                        seqs += seq[sample]
                g.write(">{}\n{}\n".format(sample, seqs))

    def get_gene_id(self, file):
        uniprot_iterator = SeqIO.parse(file, "fasta")
        records ={}
        for recordin in uniprot_iterator:
            records[recordin.id] = recordin.seq
        return records

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', metavar='[gene type]', required=True, help='nul or port')
    parser.add_argument('-s', metavar='[sample_list]', required=True, help='sample list')
    parser.add_argument('-o', metavar='[out file]', required=True, help='result fasta')
    args = parser.parse_args()
    type = args.t
    sample = args.s
    out = args.o
    ortholog = ortholog()
    ortholog.get_merge_fasta(type, sample, out)
