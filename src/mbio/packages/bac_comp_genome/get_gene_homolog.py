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
    def get_merge_fasta(self, type, homolog, clusters,  dir, sample_list, out):
        list1 = sample_list.split(",")
        cluster_list = clusters.split(",")
        with open(out, "w") as g:
            for i in list1:
                list2 = []
                if type in ['nul']:
                    list2 = self.get_gene_id(dir + "/" + i + "/" + i + "_CDS.fna")
                elif type in ['port']:
                    list2 = self.get_gene_id(dir + "/" + i + "/" + i + "_CDS.faa")
                table = pd.read_table(homolog, sep='\t', header=0)
                table = table[table['Cluster_ID'].isin(cluster_list)]
                list3 = list(table[i])
                print len(list3)
                seqs = ''
                for j in list3:
                    if re.search(",",j):
                        gene_id = j.split(",")[0].split("|")[1]
                    else:
                        gene_id =j.split("|")[1]
                    for t in list2:
                        if gene_id == t.id:
                            seqs += t.seq
                g.write(">{}\n{}\n".format(i, seqs))

    def get_gene_id(self, file):
        uniprot_iterator = SeqIO.parse(file, "fasta")
        records = list(uniprot_iterator)
        return records

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', metavar='[gene type]', required=True, help='nul or port')
    parser.add_argument('-i', metavar='[homolog cluster]', required=True, help='homolog cluster file')
    parser.add_argument('-dir', metavar='[gene dir ]', required=True, help='gene dir')
    parser.add_argument('-s', metavar='[sample_list]', required=True, help='sample list')
    parser.add_argument('-list', metavar='[cluster list]', required=True, help='clusters list')
    parser.add_argument('-o', metavar='[out file]', required=True, help='result fasta')
    args = parser.parse_args()
    type = args.t
    table = args.i
    path = args.dir
    sample = args.s
    cluster = args.list
    out = args.o
    ortholog = ortholog()
    ortholog.get_merge_fasta(type, table, cluster, path, sample, out)
