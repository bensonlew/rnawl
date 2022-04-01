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
    def get_16s_fasta(self, blast, custom, out):
        s16 ={}
        genome_fna ={}
        with open(custom, "r") as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                print(lin)
                s16[lin[0]] = lin[4]
                genome_fna[lin[0]] = lin[1]
        data = pd.read_table(blast, sep='\t', header=None)
        data.sort_values(2, inplace=True, ascending=False)
        samples =[]
        list1 =[]
        for tup in data.itertuples():
            sample =tup[2].split("_rRNA")[0]
            if sample not in samples:
                list1.append(tup[0])
                samples.append(sample)
        data = data.ix[list1,:]
        data2 = data[data[2] >= 98.7]
        sample1 = list(set(list(data[0])))[0]
        if data2.shape[0] == 0:
            with open(out + "/"+ sample1 +".taxon.xls", "w") as f:
                f.write("{}\t{}\n".format(sample1, "new species"))
        elif data2.shape[0] > 0:
            data3 = data[data[2] >= 98.7]
            data3[1], data3["bb"] = data3[1].str.split("_rRNA").str
            del data3['bb']
            data3.to_csv(out + "/all.blast_ani.xls", sep='\t', header=None, index=False)
        if data.shape[0] >= 20:
            top20 = list(data.ix[0:19, :][1])
        else:
            top20 = list(data[1])
        with open(out+"/all.16s.fa", "w") as g:
            for dd in top20:
                sample = dd.split("_rRNA")[0]
                for seq_record in SeqIO.parse(s16[sample], "fasta"):
                    if seq_record.id == dd:
                        g.write(">{}\n{}\n".format(sample, seq_record.seq))
        if os.path.exists(out+"/genome"):
            shutil.rmtree(out+"/genome")
        os.mkdir(out+"/genome")
        for dd in top20:
            sample = dd.split("_rRNA")[0]
            os.link(genome_fna[sample], out+"/genome/"+ sample+".fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', metavar='[blast file]', required=True, help='blast file')
    parser.add_argument('-r', metavar='[custom file ]', required=True, help='custom file')
    parser.add_argument('-o', metavar='[dir]', required=True, help='output dir')
    args = parser.parse_args()
    seqs = get_seq()
    seqs.get_16s_fasta(args.b, args.r, args.o)