# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os,re
import pandas as pd
import argparse
from Bio import SeqIO


class gff(object):
    """
    根据老师上传的gff文件，转化成我们需要的gff文件格式，并统计CDS，rRNA，tRNA的数量
    """
    def check_gff(self, ids, gff, out, dict_len):
        with open (gff, "r") as f, open(out + ".gff","w") as g, open(out + ".stat.xls", "w") as h:
            lines = f.readlines()
            g.write("##gff-version 3\n")
            cds_num = 0
            trna_num = 0
            rrna_num = 0
            ends =[]
            names = []
            print(ids)
            for line in lines:
                if re.search("^#", line):
                    continue
                else:
                    seq_id = ''
                    start = ''
                    end = ''
                    strand = ''
                    lin = line.strip().split("\t")
                    if lin[2] in ['region']:
                        continue
                    for i in range(0, len(lin)):
                        if lin[i] in ids:
                            print(lin[i])
                            names.append(lin[i])
                            seq_id = lin[i]
                        elif lin[i] in ["+", "-"]:
                            strand = lin[i]
                        if re.search('^[1-9]\d*$', lin[i]) and start == '':
                            start = lin[i]
                        elif re.search('^[1-9]\d*$', lin[i]) and start != '' and end == '':
                            ends.append(int(lin[i]))
                            end = lin[i]
                    if lin[1] == "RefSeq":
                        if lin[2] in ["tRNA","rRNA","CDS"]:
                            gene_id = lin[-1].split("locus_tag=")[1].split(";")[0]
                        if re.search("tRNA", lin[2]):
                            trna_num += 1
                            g.write(gene_id + "\t"+ seq_id + '\ttRNA\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" +lin[-1] + "\n")
                        elif re.search("rRNA", lin[2]):
                            rrna_num += 1
                            g.write(gene_id + "\t"+ seq_id + '\trRNA\t' + start + "\t"+end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                        elif re.search("CDS", lin[2]):
                            cds_num += 1
                            g.write(gene_id + "\t"+ seq_id + '\tCDS\t' + start + "\t"+end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                    if lin[1] in ["Genbank","Genebank"]:
                        if lin[2] in ["tRNA","rRNA","CDS"]:
                            gene_id = lin[-1].split("locus_tag=")[1].split(";")[0]
                        if re.search("tRNA", lin[2]):
                            trna_num += 1
                            g.write(gene_id + "\t"+ seq_id + '\ttRNA\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" +lin[-1] + "\n")
                        elif re.search("rRNA", lin[2]):
                            rrna_num += 1
                            g.write(gene_id + "\t"+ seq_id + '\trRNA\t' + start + "\t"+end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                        elif re.search("CDS", lin[2]):
                            cds_num += 1
                            g.write(gene_id + "\t"+ seq_id + '\tCDS\t' + start + "\t"+end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                    elif lin[1] == ".":
                        gene_id = lin[-1].split(";")[0].split("=")[1]
                        if strand == '':
                            if start > end:
                                strand = "-"
                            else:
                                strand = "+"
                        if lin[2] == "CDS":
                            cds_num += 1
                            if seq_id == '':
                                g.write(gene_id + "\t" + ids[0] + '\tCDS\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                            else:
                                g.write(gene_id + "\t" +seq_id + '\tCDS\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                        elif lin[2] == "tRNA":
                            trna_num += 1
                            if seq_id == '':
                                g.write(gene_id + "\t" + ids[0] + '\ttRNA\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                            else:
                                g.write(gene_id + "\t" +seq_id + '\ttRNA\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                        elif lin[2] == "rRNA":
                            rrna_num += 1
                            if seq_id == '':
                                g.write(gene_id + "\t" + ids[0] + '\trRNA\t' + start + "\t" + end + "\t" + "." +"\t" + strand + "\t" + "." + "\t" + lin[-1] + "\n")
                            else:
                                g.write(gene_id + "\t" +seq_id + '\trRNA\t' + start + "\t" + end + "\t" + "." +"\t"+ strand + "\t" + "." +"\t" + lin[-1] + "\n")
                    elif (lin[1] in ids) and (lin[2] in ["CDS", "rRNA", "tRNA"]) and re.search('^[1-9]\d*$', lin[3]) and re.search('^[1-9]\d*$', lin[4]) and lin[5] == "." and lin[6] in ["+", "-"]:
                        if lin[2] == "CDS":
                            cds_num += 1
                        elif lin[2] == "tRNA":
                            trna_num += 1
                        elif lin[2] == "rRNA":
                            rrna_num += 1
                        g.write("\t".join(lin)+"\n")
            ends = sorted(ends)
            names = list(set(names))
            h.write("Name\tCDS_num\ttRNA_num\trRNA_num\tgff_status\n")
            name = os.path.basename(out + ".gff")
            print(names)

            if len(names) <= 1:
                print(ends[-1])
                print(names[0])
                print(dict_len[names[0]])
                if ends[-1] > dict_len[names[0]]:
                    h.write(name + "\t"+ str(cds_num) + "\t" + str(trna_num) + "\t" + str(rrna_num) + "\t"+ "false\n")
                else:
                    h.write(name + "\t" + str(cds_num) + "\t" + str(trna_num) + "\t" + str(rrna_num) + "\t" + "true\n")
            else:
                h.write(name + "\t" + str(cds_num) + "\t" + str(trna_num) + "\t" + str(rrna_num) + "\t" + "true\n")

    def get_seqid(self, fasta):
        """
        获取序列的ids
        :param fasta:
        :return:
        """
        list =[]
        dict = {}
        iterators = SeqIO.parse(fasta, "fasta")
        for iterator in iterators:
            list.append(iterator.id)
            dict[iterator.id] = len(iterator.seq)
        return list,dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', metavar='[gff file]', required=True, help='gff file')
    parser.add_argument('-s', metavar='[genome fasta file]', required=True, help='genome fasta file')
    parser.add_argument('-o', metavar='[out file]', required=True, help='result fasta')
    args = parser.parse_args()
    type = args.g
    fa = args.s
    out = args.o
    gff = gff()
    seqids, seq_len = gff.get_seqid(fa)
    gff.check_gff(seqids, type, out, seq_len)
