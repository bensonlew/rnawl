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
    def get_merge_fasta(self, type, homolog, dir, sample_list, out):
        list1 = []
        list1 = sample_list.split(",")
        num =len(list1)
        table = pd.read_table(homolog, sep='\t', header=0)
        table = table[table['Sample_number'].isin([num]) & table['Gene_number'].isin([num])]
        rows = table.shape[0]
        dict = {}
        for i in list1:
            list2 =[]
            if type in ['nul']:
                list2 = self.get_gene_id(dir + "/" + i + "/" + i + "_CDS.fna")
            elif type in ['port']:
                list2 = self.get_gene_id(dir + "/" + i + "/" + i + "_CDS.faa")
            dict[i] = list2
        n = rows//50
        print (rows)
        if n >= 1:
            for i in range(0, n):
                table2 = table.iloc[i * 50:50 * (i + 1) - 1, :]
                print (i)
                with open(out + "/Clusters" + str(i) + ".fasta", "w") as g:
                    for j in list1:
                        seqs = ''
                        list2 = list(table2[j])
                        for w in list2:
                            geneid = w.split("|")[1]
                            if geneid in dict[j]:
                                seqs += dict[j][geneid]
                        g.write(">{}\n{}\n".format(j, seqs))
        m =rows % 50
        print ("aaaa" +str(m))
        if m != 0:
            table2 = table.iloc[n * 50:rows, :]
            with open(out + "/Clusters" + str(n) + ".fasta", "w") as g:
                for j in list1:
                    seqs = ''
                    list2 = list(table2[j])
                    for w in list2:
                        geneid = w.split("|")[1]
                        if geneid in dict[j]:
                            seqs += dict[j][geneid]
                    g.write(">{}\n{}\n".format(j, seqs))

    def get_gene_id(self, file):
        uniprot_iterator = SeqIO.parse(file, "fasta")
        records ={}
        for recordin in uniprot_iterator:
            records[recordin.id] = recordin.seq
        return records

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', metavar='[gene type]', required=True, help='nul or port')
    parser.add_argument('-i', metavar='[homolog cluster]', required=True, help='homolog cluster file')
    parser.add_argument('-dir', metavar='[gene dir ]', required=True, help='gene dir')
    parser.add_argument('-s', metavar='[sample_list]', required=True, help='sample list')
    parser.add_argument('-o', metavar='[out file]', required=True, help='result fasta')
    args = parser.parse_args()
    type = args.t
    table = args.i
    path = args.dir
    sample = args.s
    out = args.o
    ortholog = ortholog()
    ortholog.get_merge_fasta(type, table, path, sample, out)
