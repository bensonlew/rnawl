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
    根据blast的比对结果提取样品的16s序列以及基因组文件
    """
    def get_16s_fasta(self, blast, s16, out):
        data = pd.read_table(blast, sep='\t', header=None)
        data.sort_values(2, inplace=True, ascending=False)
        if data.shape[0] >= 20:
            top20 = list(data.ix[0:19, :][1])
        else:
            top20 = list(data[1])
        with open(out, "w") as g:
            for seq_record in SeqIO.parse(s16, "fasta"):
                if seq_record.id in top20:
                    name = "_".join(seq_record.id.split("_")[1:])
                    g.write(">{}\n{}\n".format(name, seq_record.seq))

    def get_genome_fasta(self, blast, dir, out):
        if os.path.exists(out+"/genome"):
            shutil.rmtree(out+"/genome")
        os.mkdir(out+"/genome")
        data = pd.read_table(blast, sep='\t', header=None)
        data.sort_values(2, inplace=True, ascending=False)
        if data.shape[0] >=20:
            top20 = list(data.ix[0:19, :][1])
        else:
            top20 = list(data[1])
        data2 = data[data[2] >= 98.7]
        sample1 = list(set(list(data[0])))[0]
        if data2.shape[0] == 0:
            with open(out + "/"+ sample1 +".taxon.xls", "w") as f:
                f.write("{}\t{}\n".format(sample1, "new species"))
        elif data2.shape[0] > 0:
            data2.to_csv(out + "/all.blast_ani.xls", sep='\t', header=None, index=False)
        for dd in top20:
            sample ="_".join(dd.split("_")[1:])
            path = dir +"/"+sample+"_genomic.fna.gz"
            if os.path.exists(path):
                fna = out+"/genome/" + sample + ".fasta"
                f = gzip.open(path, "rb")
                f_content = f.read()
                f.close()
                open(fna, "wb+").write(f_content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', metavar='[blast file]', required=True, help='blast file')
    parser.add_argument('-r', metavar='[rrna file ]', required=True, help='16s rRNA fasta file')
    parser.add_argument('-d', metavar='[genomes dir]', required=True, help='genomes dir')
    parser.add_argument('-o', metavar='[dir]', required=True, help='output dir')
    args = parser.parse_args()
    seqs = get_seq()
    seqs.get_16s_fasta(args.b, args.r, args.o+"/all.16s.fa")
    seqs.get_genome_fasta(args.b, args.d, args.o)