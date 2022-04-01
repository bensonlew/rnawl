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
    def get_genome_fasta(self, blast, ncbi, out):
        if os.path.exists(out+"/genome"):
            shutil.rmtree(out+"/genome")
        os.mkdir(out+"/genome")
        data = pd.read_table(blast, sep='\t', header=None)
        top20 = list(data[1])
        with open(ncbi, "r") as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                if lin[0] in top20:
                    os.link(lin[1], out+"/genome/"+ lin[0]+".fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', metavar='[blast file]', required=True, help='blast file')
    parser.add_argument('-d', metavar='[download NCBI file]', required=True, help='download NCBI file')
    parser.add_argument('-o', metavar='[dir]', required=True, help='output dir')
    args = parser.parse_args()
    seqs = get_seq()
    seqs.get_genome_fasta(args.b, args.d, args.o)