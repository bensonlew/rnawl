# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os,re
import pandas as pd
import argparse
from Bio import SeqIO


class ortholog(object):
    """
    antismash的丰度统计和gene_list生成
    """
    def get_merge_fasta(self, coregene, coregenes, sample_list, out):
        list1 = sample_list.split(",")
        coregene_list = coregenes.split(",")
        with open(out, "w") as g:
            for i in list1:
                table = pd.read_table(coregene, sep='\t', header=0)
                table = table[table['Bin_id'].isin([i]) & table['Name'].isin(coregene_list)]
                list3 = list(table['Protein'])
                seqs = ''
                print len(list3)
                if len(list3):
                    for j in list3:
                        seqs += j
                    g.write(">{}\n{}\n".format(i, seqs))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[core_gene blast]', required=True, help='core_gene blast file')
    parser.add_argument('-s', metavar='[sample_list]', required=True, help='sample list')
    parser.add_argument('-list', metavar='[houskeeping name list]', required=True, help='houskeeping name list')
    parser.add_argument('-o', metavar='[out file]', required=True, help='result fasta')
    args = parser.parse_args()
    table = args.i
    sample = args.s
    cluster = args.list
    out = args.o
    ortholog = ortholog()
    ortholog.get_merge_fasta(table, cluster, sample, out)
