# -*- coding: utf-8 -*-
# __author__ : 'shicaiping'
# __date__: 20200305

import pandas as pd
from Bio import SeqIO
import os


class ExpUnitsConvert(object):
    """
    Used for convert RNA-seq expression units convert, count/CPM/RPM/TPM/FPKM/RPKM are implemented.

    RPM/CPM: Reads/Counts of exon model per million mapped reads
    RPM/CPM=Total exon reads/ Mapped reads(Millions)

    RPKM/FPKM: Reads/Fragments Per Kilobase of exon model per Million mapped reads
    RPKM/FPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)]

    TPM is like RPKM and FPKM, except the order of operation is switched.
    """

    def __init__(self, exp_matrix, gene_length, gene_fasta, convert_type, intersect=True, num=2):
        """
        :param count: Sample name is used as column name, first column must be gene id
        :param tpm: Sample name is used as column name, first column must be gene id
        :param cpm: Sample name is used as column name, first column must be gene id
        :param fpkm: Sample name is used as column name, first column must be gene id
        :param gene_length: This file must contains at least two columns, First column must be gene id, second column must be gene length.
        :param convert: count2tpm, count2cpm, count2fpkm, fpkm2tmp, cpm2tpm, cpm2fpkm, fpkm2cpm
        :return: converted gene expression matrix
        """
        convert_type_legal = ["count2tpm", "count2cpm", "count2fpkm", "fpkm2tpm", "cpm2tpm", "cpm2fpkm", "fpkm2cpm"]
        if not convert_type:
            raise Exception("Please input convert type parameter")
        if convert_type not in convert_type_legal:
            raise Exception("This convert type cannot be supported!")
        if convert_type in ["count2tpm", "count2cpm", "count2fpkm"]:
            self.count = exp_matrix
        if convert_type in ["cpm2tpm", "cpm2fpkm"]:
            self.cpm = exp_matrix
        if convert_type in ["fpkm2tpm", "fpkm2cpm"]:
            self.fpkm = exp_matrix
        if gene_length:
            self.gene_length = gene_length
        elif gene_fasta:
            self.gene_length = self.fasta_length(gene_fasta)
        else:
            self.gene_length = None
        self.convert = convert_type
        self.num = num
        self.intersect = intersect

    def run(self):
        if self.convert == "count2cpm":
            if not self.count:
                raise Exception("You must be input count matrix to finish {}!".format(self.convert))
            self.count2cpm(count=self.count, num=self.num)
        if self.convert == "count2tpm":
            if not self.count:
                raise Exception("You must be input count matrix to finish {}!".format(self.convert))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix  or gene fasta file to finish {}!".format(self.convert))
            self.count2tpm(count=self.count, gene_length=self.gene_length, intersect=self.intersect, num=self.num)
        if self.convert == "count2fpkm":
            if not self.count:
                raise Exception("You must be input count matrix to finish {}!".format(self.convert))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix  or gene fasta file to finish {}!".format(self.convert))
            self.count2fpkm(count=self.count, gene_length=self.gene_length, intersect=self.intersect, num=self.num)
        if self.convert == "fpkm2tpm":
            if not self.fpkm:
                raise Exception("You must be input fpkm/rpkm matrix to finish {}!".format(self.convert))
            self.fpkm2tpm(fpkm=self.fpkm, num=self.num)
        if self.convert == "cpm2fpkm":
            if not self.cpm:
                raise Exception("You must be input cpm matrix to finish {}!".format(self.convert))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix  or gene fasta file to finish {}!".format(self.convert))
            self.cpm2fpkm(cpm=self.cpm, gene_length=self.gene_length, intersect=self.intersect, num=self.num)
        if self.convert == "cpm2tpm":
            if not self.cpm:
                raise Exception("You must be input cpm matrix to finish {}!".format(self.convert))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix  or gene fasta file to finish {}!".format(self.convert))
            self.cpm2tpm(cpm=self.cpm, gene_length=self.gene_length, intersect=self.intersect, num=self.num)
        if self.convert == "fpkm2cpm":
            if not self.fpkm:
                raise Exception("You must be input fpkm/rpkm matrix to finish {}!".format(self.convert))
            if not self.gene_length:
                raise Exception(
                    "You must be input gene length matrix or gene fasta file to finish {}!".format(self.convert))
            self.fpkm2cpm(fpkm=self.fpkm, gene_length=self.gene_length, intersect=self.intersect, num=self.num)

    def fasta_length(self, fasta):
        length_file = "gene_length.txt"
        with open(length_file, "w") as w:
            w.write("Gene_ID\tLength\n")
            for seq_record in SeqIO.parse(fasta, "fasta"):
                gene_id = seq_record.id
                gene_length = len(seq_record.seq)
                w.write(gene_id + "\t" + str(gene_length) + "\n")
        return length_file

    def count2cpm(self, count, num):
        count_table = pd.read_table(count, index_col=0)
        columns = count_table.columns
        cpm_dict = dict()
        for sample in columns:
            total_counts = sum(count_table[sample])
            cpm_dict[sample] = (count_table[sample] / total_counts * 1000000).round(num)
        df_cpm = pd.DataFrame(cpm_dict)
        df_cpm.to_csv(os.path.basename(count) + '.count2cpm.xls', sep='\t')

    def count2tpm(self, count, gene_length, intersect, num):
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
        df_tpm.to_csv(os.path.basename(count) + '.count2tpm.xls', sep='\t')

    def count2fpkm(self, count, gene_length, intersect, num):
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
        df_fpkm.to_csv(os.path.basename(count) + '.count2fpkm.xls', sep='\t')

    def fpkm2tpm(self, fpkm, num):
        fpkm_table = pd.read_table(fpkm, index_col=0)
        tpm_dict = dict()
        columns = fpkm_table.columns
        for sample in columns:
            total_fpkm = sum(fpkm_table[sample])
            tpm = (fpkm_table[sample] / total_fpkm * 1000000).round(num)
            tpm_dict[sample] = tpm
        df_tpm = pd.DataFrame(tpm_dict)
        df_tpm.to_csv(os.path.basename(fpkm) + '.fpkm2tpm.xls', sep='\t')

    def cpm2fpkm(self, cpm, gene_length, intersect, num):
        cpm_table = pd.read_table(cpm, index_col=0)
        gene_length = pd.read_table(gene_length, index_col=0)
        ## 去除基因长度为0的行
        gene_length = gene_length[gene_length[gene_length.columns[0]] > 0]
        ## 判断gene id是否一致，如果不一致，取交集或报错
        if set(gene_length.index).difference(set(cpm_table.index)) or set(cpm_table.index).difference(
                set(gene_length.index)):
            if intersect:
                print "cpm table and gene length table have different gene ids, choose the intersection part."
                intersect_ids = set(gene_length.index).intersection(set(cpm_table.index))
                print "there are {} intersect gene ids.".format(len(intersect_ids))
                if len(intersect_ids) == 0:
                    raise Exception("There are no intersect gene ids.")
                gene_length = gene_length.loc[intersect_ids]
                cpm_table = cpm_table.loc[intersect_ids]
            else:
                raise Exception("cpm table and gene length table have different gene ids, such as {}, {}".format(
                    set(gene_length.index).difference(set(cpm_table.index)),
                    set(cpm_table.index).difference(set(gene_length.index))))
        cpm_table.sort_index(inplace=True)
        gene_length.sort_index(inplace=True)
        columns = cpm_table.columns
        gene_len = gene_length[gene_length.columns[0]]
        fpkm_dict = dict()
        for sample in columns:
            fpkm = (cpm_table[sample] / gene_len * 1000).round(num)
            fpkm_dict[sample] = fpkm
        df_fpkm = pd.DataFrame(fpkm_dict)
        df_fpkm.to_csv(os.path.basename(cpm) + '.cpm2fpkm.xls', sep='\t')

    def cpm2tpm(self, cpm, gene_length, intersect, num):
        cpm_table = pd.read_table(cpm, index_col=0)
        gene_length = pd.read_table(gene_length, index_col=0)
        ## 去除基因长度为0的行
        gene_length = gene_length[gene_length[gene_length.columns[0]] > 0]
        ## 判断gene id是否一致，如果不一致，取交集或报错
        if set(gene_length.index).difference(set(cpm_table.index)) or set(cpm_table.index).difference(
                set(gene_length.index)):
            if intersect:
                print "cpm table and gene length table have different gene ids, choose the intersection part."
                intersect_ids = set(gene_length.index).intersection(set(cpm_table.index))
                print "there are {} intersect gene ids.".format(len(intersect_ids))
                if len(intersect_ids) == 0:
                    raise Exception("There are no intersect gene ids.")
                gene_length = gene_length.loc[intersect_ids]
                cpm_table = cpm_table.loc[intersect_ids]
            else:
                raise Exception("cpm table and gene length table have different gene ids, such as {}, {}".format(
                    set(gene_length.index).difference(set(cpm_table.index)),
                    set(cpm_table.index).difference(set(gene_length.index))))
        cpm_table.sort_index(inplace=True)
        gene_length.sort_index(inplace=True)
        columns = cpm_table.columns
        gene_len = gene_length[gene_length.columns[0]]
        fpkm_dict = dict()
        for sample in columns:
            fpkm = (cpm_table[sample] / gene_len * 1000).round(num)
            fpkm_dict[sample] = fpkm
        df_fpkm = pd.DataFrame(fpkm_dict)
        tpm_dict = dict()
        columns = df_fpkm.columns
        for sample in columns:
            total_fpkm = sum(df_fpkm[sample])
            tpm = (df_fpkm[sample] / total_fpkm * 1000000).round(num)
            tpm_dict[sample] = tpm
        df_tpm = pd.DataFrame(tpm_dict)
        df_tpm.to_csv(os.path.basename(cpm) + '.cpm2tpm.xls', sep='\t')

    def fpkm2cpm(self, fpkm, gene_length, intersect, num):
        fpkm_table = pd.read_table(fpkm, index_col=0)
        gene_length = pd.read_table(gene_length, index_col=0)
        ## 去除基因长度为0的行
        gene_length = gene_length[gene_length[gene_length.columns[0]] > 0]
        ## 判断gene id是否一致，如果不一致，取交集或报错
        if set(gene_length.index).difference(set(fpkm_table.index)) or set(fpkm_table.index).difference(
                set(gene_length.index)):
            if intersect:
                print "fpkm table and gene length table have different gene ids, choose the intersection part."
                intersect_ids = set(gene_length.index).intersection(set(fpkm_table.index))
                print "there are {} intersect gene ids.".format(len(intersect_ids))
                if len(intersect_ids) == 0:
                    raise Exception("There are no intersect gene ids.")
                gene_length = gene_length.loc[intersect_ids]
                fpkm_table = fpkm_table.loc[intersect_ids]
            else:
                raise Exception("fpkm table and gene length table have different gene ids, such as {}, {}".format(
                    set(gene_length.index).difference(set(fpkm_table.index)),
                    set(fpkm_table.index).difference(set(gene_length.index))))
        fpkm_table.sort_index(inplace=True)
        gene_length.sort_index(inplace=True)
        columns = fpkm_table.columns
        gene_len = gene_length[gene_length.columns[0]]
        cpm_dict = dict()
        for sample in columns:
            cpm = (fpkm_table[sample] * gene_len / 1000).round(num)
            cpm_dict[sample] = cpm
        df_cpm = pd.DataFrame(cpm_dict)
        df_cpm.to_csv(os.path.basename(fpkm) + '.fpkm2cpm.xls', sep='\t')


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='This script is used to convert rna-seq expression units.')
    parser.add_argument('-convert_type', type=str,
                        help='exp units to be convert, can be one of these: count2tpm, count2cpm, count2fpkm, fpkm2tmp, cpm2tpm, cpm2fpkm, fpkm2cpm')
    parser.add_argument('-exp_matrix', type=str,
                        help='path of exp matrix table. Sample name is used as column name, first column must be gene id')
    parser.add_argument('-gene_length', type=str,
                        help='path of gene length table. This file must contains at least two columns, First column must be gene id, second column must be gene length.')
    parser.add_argument('-gene_fasta', type=str,
                        help='path of gene fasta table. When gene length file is not provided, this file can be used to compute gene length.')
    parser.add_argument('-num', type=int, default=2, help='Number of decimal places reserved.')
    parser.add_argument('-intersect', type=bool, default=True,
                        help='when gene ids are not complete the same, choose the intersection part.')

    args = parser.parse_args()

    convert = ExpUnitsConvert(args.exp_matrix, args.gene_length, args.gene_fasta,
                              args.convert_type, intersect=args.intersect, num=args.num)
    convert.run()
