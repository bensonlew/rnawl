# -*- coding: utf-8 -*-
# __author__ : 'shicaiping,zoujiaxun'
# __date__: 20200426

import pandas as pd
from Bio import SeqIO
import os
import subprocess
import logging

class ExpNorm(object):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.

    TMM(Trimmed Mean of M) is a normalized method for edgeR
    DEGseq2 is a normalized method for DEGseq2
    uqua is upperquartile
    """

    def __init__(self, exp_matrix, gene_length, gene_fasta, method, output, intersect=True, num=2):
        """
        :param count: Sample name is used as column name, first column must be gene id
        :param tpm: Sample name is used as column name, first column must be gene id
        :param cpm: Sample name is used as column name, first column must be gene id
        :param fpkm: Sample name is used as column name, first column must be gene id
        :param gene_length: This file must contains at least two columns, First column must be gene id, second column must be gene length.
        :param convert: count2tpm, count2cpm, count2fpkm, fpkm2tmp, cpm2tpm, cpm2fpkm, fpkm2cpm
        :return: converted gene expression matrix
        """
        method_legal = ["tpm", "cpm", "fpkm", 'TMM', 'TMMwzp', 'RLF', 'uqua', 'DESeq2']
        if not method:
            raise Exception("Please input method parameter")
        if method not in method_legal:
            raise Exception("This convert type cannot be supported!")

        if gene_length:
            self.gene_length = gene_length
        elif gene_fasta:
            self.gene_length = self.fasta_length(gene_fasta)
        else:
            self.gene_length = None
        self.count = exp_matrix
        self.method = method
        self.num = num
        self.intersect = intersect
        self.output = output

    def run(self):
        if self.method == "cpm":
            if not self.count:
                raise Exception("You must be input count matrix to finish {}!".format(self.method))
            self.cpm(count=self.count, num=self.num)
        if self.method == "tpm":
            if not self.count:
                raise Exception("You must be input count matrix to finish {}!".format(self.method))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix  or gene fasta file to finish {}!".format(self.method))
            self.tpm(count=self.count, gene_length=self.gene_length, intersect=self.intersect, num=self.num)
        if self.method == "fpkm":
            if not self.count:
                raise Exception("You must be input count matrix to finish {}!".format(self.method))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix  or gene fasta file to finish {}!".format(self.method))
            self.fpkm(count=self.count, gene_length=self.gene_length, intersect=self.intersect, num=self.num)
        if self.method in ['TMM', 'TMMwzp', 'RLE', 'uqua', 'DESeq2']:
            cmd = '{} {} -i {} -m {} -n {} -o {}'.format(
                args.r_interpreter, args.r_script, self.count, self.method, self.num, args.output)
            proc = subprocess.Popen(
                cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            retcode = proc.wait()
            outs, errs = proc.communicate()
            if retcode:
                msg = '\n'.join(
                    ('fail to excecute {}'.format(cmd), 'STDOUT: {}'.format(outs), 'STDERR: {}'.format(errs)))
                raise Exception(msg)
            else:
                logging.info('succeed in excecuting command ({})'.format(cmd))
                logging.debug(outs)
                logging.error(errs)



    def fasta_length(self, fasta):
        length_file = "gene_length.txt"
        with open(length_file, "w") as w:
            w.write("Gene_ID\tLength\n")
            for seq_record in SeqIO.parse(fasta, "fasta"):
                gene_id = seq_record.id
                gene_length = len(seq_record.seq)
                w.write(gene_id + "\t" + str(gene_length) + "\n")
        return length_file

    def cpm(self, count, num):
        count_table = pd.read_table(count, index_col=0)
        columns = count_table.columns
        cpm_dict = dict()
        for sample in columns:
            total_counts = sum(count_table[sample])
            cpm_dict[sample] = (count_table[sample] / total_counts * 1000000).round(num)
        df_cpm = pd.DataFrame(cpm_dict)
        df_cpm.to_csv(self.output, sep='\t')

    def tpm(self, count, gene_length, intersect, num):
        count_table = pd.read_table(count, index_col=0)
        gene_length = pd.read_table(gene_length, index_col=0)
        ## 去除基因长度为0的行
        gene_length = gene_length[gene_length[gene_length.columns[0]] > 0]
        ## 判断gene id是否一致，如果不一致，取交集或报错
        if set(gene_length.index).difference(set(count_table.index)) or set(count_table.index).difference(
                set(gene_length.index)):
            if intersect:
                print "count table and gene length table have different gene ids, choose the intersection part."
                intersect_ids = set(gene_length.index).intersection(set(count_table.index))
                print "there are {} intersect gene ids.".format(len(intersect_ids))
                if len(intersect_ids) == 0:
                    raise Exception("There are no intersect gene ids.")
                gene_length = gene_length.loc[intersect_ids]
                count_table = count_table.loc[intersect_ids]
            else:
                raise Exception("count table and gene length table have different gene ids, such as {}, {}".format(
                    set(gene_length.index).difference(set(count_table.index)),
                    set(count_table.index).difference(set(gene_length.index))))
        count_table.sort_index(inplace=True)
        gene_length.sort_index(inplace=True)
        columns = count_table.columns
        gene_len = gene_length[gene_length.columns[0]]
        tpm_dict = dict()
        for sample in columns:
            rpk = count_table[sample] / gene_len
            norm_gene_len_total_counts = sum(rpk)
            tpm = (rpk / norm_gene_len_total_counts * 1000000).round(num)
            tpm_dict[sample] = tpm
        df_tpm = pd.DataFrame(tpm_dict)
        df_tpm.to_csv(self.output, sep='\t')

    def fpkm(self, count, gene_length, intersect, num):
        count_table = pd.read_table(count, index_col=0)
        gene_length = pd.read_table(gene_length, index_col=0)
        ## 去除基因长度为0的行
        gene_length = gene_length[gene_length[gene_length.columns[0]] > 0]
        ## 判断gene id是否一致，如果不一致，取交集或报错
        if set(gene_length.index).difference(set(count_table.index)) or set(count_table.index).difference(
                set(gene_length.index)):
            if intersect:
                print "count table and gene length table have different gene ids, choose the intersection part."
                intersect_ids = set(gene_length.index).intersection(set(count_table.index))
                print "there are {} intersect gene ids.".format(len(intersect_ids))
                if len(intersect_ids) == 0:
                    raise Exception("There are no intersect gene ids.")
                gene_length = gene_length.loc[intersect_ids]
                count_table = count_table.loc[intersect_ids]
            else:
                raise Exception("count table and gene length table have different gene ids, such as {}, {}".format(
                    set(gene_length.index).difference(set(count_table.index)),
                    set(count_table.index).difference(set(gene_length.index))))
        count_table.sort_index(inplace=True)
        gene_length.sort_index(inplace=True)
        columns = count_table.columns
        gene_len = gene_length[gene_length.columns[0]]
        fpkm_dict = dict()
        for sample in columns:
            rpk = count_table[sample] / gene_len
            total_counts = sum(count_table[sample])
            fpkm = (rpk / total_counts * 1000000 * 1000).round(num)
            fpkm_dict[sample] = fpkm
        df_fpkm = pd.DataFrame(fpkm_dict)
        df_fpkm.to_csv(self.output, sep='\t')


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='This script is used to convert rna-seq expression units.')
    parser.add_argument('-convert_type', type=str,
                        help='exp norm_method, can be one of these: tpm, cpm, fpkm')
    parser.add_argument('-exp_matrix', type=str,
                        help='path of exp matrix table. Sample name is used as column name, first column must be gene id')
    parser.add_argument('-gene_length', type=str,
                        help='path of gene length table. This file must contains at least two columns, First column must be gene id, second column must be gene length.')
    parser.add_argument('-gene_fasta', type=str,
                        help='path of gene fasta table. When gene length file is not provided, this file can be used to compute gene length.')
    parser.add_argument('-num', type=int, default=2, help='Number of decimal places reserved.')
    parser.add_argument('-intersect', type=bool, default=True,
                        help='when gene ids are not complete the same, choose the intersection part.')
    parser.add_argument('-output', type=str,
                        help='path of exp_norm table')
    parser.add_argument('--r_interpreter', action='store',
                          help='path of R interpreter', metavar='<FILE>', dest='r_interpreter')
    parser.add_argument('--r_script', action='store',
                          help='requried R script', metavar='<FILE>', dest='r_script')

    args = parser.parse_args()

    convert = ExpNorm(args.exp_matrix, args.gene_length, args.gene_fasta,
                              args.convert_type, intersect=args.intersect, num=args.num, output=args.output)
    convert.run()
