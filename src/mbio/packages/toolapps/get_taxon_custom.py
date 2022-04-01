# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.06

import os,re
import pandas as pd
import argparse
from Bio import SeqIO
import gzip
import shutil

class get_seq(object):
    """
    根据blast的比对结果提取样品的ani基因组文件
    """
    def get_genome_fasta(self, query, blast, ani, database, out):
        tson ={}
        with open(database, "r") as g:
            lines = g.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                tson[lin[0]] = lin[2]
        blast1 = {}
        with open(blast, "r") as d:
            lines = d.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                blast1[lin[1]] = lin[2]
        data = pd.read_table(ani, sep=' ', header=None)
        data.sort_values(2, inplace=True, ascending=False)
        top20 = list(data[1])
        data2 = data[data[2] >= 95]
        with open(out+".taxon.xls", "w") as f,open(out+".anno.xls", "w") as g:
            g.write("sample\tref\tsimilarity\tani\tNCBI taxon\n")
            if data2.shape[0] ==0:
                f.write("{}\t{}\n".format(query, "new species"))
                for dd in top20:
                    sample = dd.split("/")[-1].split(".fasta")[0]
                    ani = data.ix[top20.index(dd),2]
                    g.write("{}\t{}\t{}\t{}\t{}\n".format(query, sample, str(blast[sample]), str(ani), tson[sample]))
            elif data2.shape[0] >0:
                name = data2.ix[0,1].split("/")[-1].split(".fasta")[0]
                f.write("{}\t{}\n".format(query, tson[name]))
                for ss in top20:
                    sample = ss.split("/")[-1].split(".fasta")[0]
                    ani = data.ix[top20.index(ss), 2]
                    g.write("{}\t{}\t{}\t{}\t{}\n".format(query, sample, str(blast1[sample]), str(ani), tson[sample]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='[sample name]', required=True, help='sample name')
    parser.add_argument('-b', metavar='[blast file]', required=True, help='blast file')
    parser.add_argument('-a', metavar='[ani file]', required=True, help='ani file')
    parser.add_argument('-t', metavar='[taxon]', required=True, help='database taxon file')
    parser.add_argument('-o', metavar='[dir]', required=True, help='output dir')
    args = parser.parse_args()
    seqs = get_seq()
    seqs.get_genome_fasta(args.s, args.b, args.a, args.t, args.o)