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
    def get_genome_fasta(self, blast, dir, out):
        if os.path.exists(out+"/genome"):
            shutil.rmtree(out+"/genome")
        os.mkdir(out+"/genome")
        data = pd.read_table(blast, sep='\t', header=None)
        top20 = list(data[1])
        for dd in top20:
            sample ="_".join(dd.split("_")[1:])
            path = dir +"/"+sample+"_genomic.fna.gz"
            if os.path.exists(path):
                fna = out+"/genome/" + sample + "_genomic.fna"
                f = gzip.open(path, "rb")
                f_content = f.read()
                f.close()
                open(fna, "wb+").write(f_content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', metavar='[blast file]', required=True, help='blast file')
    parser.add_argument('-d', metavar='[genomes dir]', required=True, help='genomes dir')
    parser.add_argument('-o', metavar='[dir]', required=True, help='output dir')
    args = parser.parse_args()
    seqs = get_seq()
    seqs.get_genome_fasta(args.b, args.d, args.o)