# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.01.06

import os,re
import pandas as pd
import argparse
from Bio import SeqIO


class get_seq(object):
    """
    通过gff获取序列
    """
    def get_merge_fasta(self, genome, gff, s16, faa, gff1):
        dict_faa = []
        dict_16s = []
        with open(gff, "r") as f, open(gff1, "w") as g:
            g.write("##gff-version 3\n")
            lines = f.readlines()
            for line in lines:
                if re.search("^#", line):
                    continue
                else:
                    lin = line.strip().split("\t")
                    if lin[2] in ['CDS', 'cds']:
                        gene_id =lin[8].split(";")[0].split("ID=")[1]
                        dict_faa.append([gene_id, lin[0], lin[3], lin[4], lin[6]])
                        g.write("{}\n".format("\t".join([gene_id, lin[0], lin[2],lin[3],lin[4],lin[5],lin[6],lin[7],lin[8]])))
                    elif lin[2] in ['rRNA', 'rrna']:
                        if re.search("16S ribosomal RNA", lin[8]):
                            dict_16s.append([gene_id, lin[0], lin[3], lin[4], lin[6]])
        if len(dict_16s) >0:
            sample =os.path.basename(s16).split(".16s.fasta")[0]
            self.get_16s(genome, s16, dict_16s, 1, 2, 3, sample)
        if len(dict_faa) >0:
            self.get_seq(genome, faa, dict_faa, 0, 1, 2, 3, 4)

    def get_16s(self, fa, out, list1, loction, start, end, sample):
        with open(out, "w") as g:
            num = 0
            for i in list1:
                num +=1
                id = sample+"_rRNA"+str(num)
                for seq_record in SeqIO.parse(fa, "fasta"):
                    if i[loction] == seq_record.id:
                        if int(i[end]) < int(i[start]):
                            seq1 = str(seq_record.seq[int(i[end]) - 1:int(i[start])])
                            seq = self.dna_reverse(self.dna_complement(seq1))
                            g.write(">{}\n{}\n".format(id, seq))
                        elif int(i[end]) > int(i[start]):
                            seq = str(seq_record.seq[int(i[start]) - 1:int(i[end])])
                            g.write(">{}\n{}\n".format(id, seq))
                    else:
                        raise ValueError("gff的location与基因组的location不一致！")

    def get_seq(self, fa, out, list1, gene_id, loction, start, end, strand):
        with open(out, "w") as g:
            for i in list1:
                for seq_record in SeqIO.parse(fa, "fasta"):
                    if i[loction] == seq_record.id:
                        if i[strand] == "+":
                            seq = str(seq_record.seq[int(i[start]) - 1:int(i[end])])
                            g.write(">{}\n{}\n".format(i[gene_id], seq))
                        elif i[strand] == "-":
                            if int(i[end]) < int(i[start]):
                                seq1 = str(seq_record.seq[int(i[end]) - 1:int(i[start])])
                                seq = self.dna_reverse(self.dna_complement(seq1))
                                g.write(">{}\n{}\n".format(i[gene_id], seq))
                            elif int(i[end]) > int(i[start]):
                                seq1 = str(seq_record.seq[int(i[start]) - 1:int(i[end])])
                                seq = self.dna_reverse(self.dna_complement(seq1))
                                g.write(">{}\n{}\n".format(i[gene_id], seq))
                        elif i[strand] == ".":
                            if int(i[end]) < int(i[start]):
                                seq1 = str(seq_record.seq[int(i[end]) - 1:int(i[start])])
                                seq = self.dna_reverse(self.dna_complement(seq1))
                                g.write(">{}\n{}\n".format(i[gene_id], seq))
                            elif int(i[end]) > int(i[start]):
                                seq = str(seq_record.seq[int(i[start]) - 1:int(i[end])])
                                g.write(">{}\n{}\n".format(i[gene_id], seq))


    def dna_complement(self, sequence):
        sequence = sequence.upper()
        sequence = sequence.replace('A', 't')
        sequence = sequence.replace('T', 'a')
        sequence = sequence.replace('C', 'g')
        sequence = sequence.replace('G', 'c')
        return sequence.upper()

    def dna_reverse(self, sequence):
        sequence = sequence.upper()
        return sequence[::-1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='[genome seq]', required=True, help='genome seq file')
    parser.add_argument('-g', metavar='[gff file ]', required=True, help='genome gff file')
    parser.add_argument('-p', metavar='[sample name]', required=True, help='prefix of output files')
    args = parser.parse_args()
    seqs = get_seq()
    seqs.get_merge_fasta(args.s, args.g, args.p+".16s.fasta", args.p+".CDS.fasta", args.p+".CDS.gff")