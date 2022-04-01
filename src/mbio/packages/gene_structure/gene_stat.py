# -*- coding: utf-8 -*-
# __author__ = 'guhaidong-20170829'

import os, sys
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='file')
parser.add_argument("-gene_dir", "--gene_dir", required=True)
parser.add_argument("-output_stat", "--output_stat", required=True)
parser.add_argument("-output_fa", "--output_fa", required=True)
parser.add_argument("-sp_name", "--sp_name", default=0, type=int)
args = vars(parser.parse_args())

gene_dir = args["gene_dir"]
stat = args["output_stat"]
fasta = args["output_fa"]
has_sp_name = args["sp_name"]
if os.path.exists(fasta):
    os.remove(fasta)
if os.path.exists(fasta + '_mix'):
    os.remove(fasta + '_mix')
if os.path.exists(fasta + '_sample'):
    os.remove(fasta + '_sample')


def get_gene_stat(i_gene_dir, o_stat, o_fasta):
    all_files = os.listdir(i_gene_dir)
    name_list = []
    name_dic = {}
    all_stat = dict()
    mode = 1
    for files in all_files:
        name_dump = os.path.basename(files).partition(".metagene")[0]
        if name_dump not in name_list:
            name_list.append(name_dump)
            name_dic[name_dump] = files
    if "Megahit_Mix" in name_list or "newbler" in name_list:
        mode = 2
    else:
        mode = 1
    for name in name_list:
        # stat_tmp = StatFile(i_gene_dir + "/" + files, o_fasta, name_dump, mode)
        stat_tmp = StatFile(i_gene_dir + "/" + name_dic[name], o_fasta, name, mode)
        if name not in all_stat:
            all_stat[name] = stat_tmp.static
    with open(o_stat, "w") as fw:
        fw.write("#Sample\tORFs\tTotal length(bp)\tAverage length(bp)\tMax(bp)\tMin(bp)\n")
        total_orfs = 0
        total_bases = 0
        orf_max = 0
        orf_min = 1000000
        for one in name_list:
            one_list = all_stat[one]
            # stat_line = "\t".join(one_list[0:-1])
            stat_line = "\t".join(one_list)
            fw.write(one + "\t" + stat_line + "\n")
            total_orfs += int(one_list[0])
            total_bases += int(one_list[1])
            if orf_max < int(one_list[3]):
                orf_max = int(one_list[3])
            if orf_min > int(one_list[4]):
                orf_min = int(one_list[4])
        orf_average = float(total_bases) / float(total_orfs)
        fw.write('Total' + '\t' + str(total_orfs) + '\t' + str(total_bases) + '\t' + str(orf_average) + '\t' +
                 str(orf_max) + '\t' + str(orf_min))  # 输入统计信息


class StatFile(object):
    def __init__(self, gene, o_fasta, name, mode):  # change param
        self.gene = gene
        self.name = name
        self.output = o_fasta
        self.mode = mode
        # self.n50 = 0
        # self.n90 = 0
        self.max = 0
        self.min = 1000000
        self.reads_num = 0
        self.base_num = 0
        self.ave = 0
        # self.length_list = []
        # self.limit = min_gene
        self.get_stat()
        # self.static = [str(self.reads_num), str(self.base_num), str(self.n50), str(self.n90), str(self.max), str(self.min), str(kmer)]
        self.static = [str(self.reads_num), str(self.base_num), str(self.ave), str(self.max), str(self.min)]

    def get_stat(self):
        """
        对类计算reads数，base数，max_len，min_len
        :return:
        """
        if self.mode == 2:
            fw3 = open(self.output, "a")
            if self.name in ["newbler", "Megahit_Mix"]:
                out = self.output + "_mix"
            else:
                out = self.output + "_sample"
        else:
            out = self.output
        seq = ""
        with open(self.gene) as fr, open(out, "a") as fw2:
            for line in fr:
                line = line.strip()
                if line.startswith(">"):
                    new_line = line.lstrip(">")
                    new_line = new_line if has_sp_name else self.name + '_' + new_line
                    if seq == "":
                        fw2.write('>' + new_line + '\n')
                        if self.mode == 2:
                            fw3.write('>' + new_line + '\n')
                    else:
                        fw2.write(seq + '\n>' + new_line + '\n')
                        if self.mode == 2:
                            fw3.write(seq + '\n>' + new_line + '\n')
                        self.calculate(seq)
                        seq = ""
                else:
                    seq += line
            fw2.write(seq + '\n')
            if self.mode == 2:
                fw3.write(seq + '\n')
                fw3.close()
            self.calculate(seq)  # 对最后一个序列进行计算
            self.ave = float(self.base_num) / float(self.reads_num)

    def calculate(self, sequence):
        """
        计算最大片段，最小片段，总长，总reads数
        :return:
        """
        self.reads_num += 1
        self.base_num += len(sequence)
        # self.length_list.append(len(sequence))
        if self.max < len(sequence):
            self.max = len(sequence)
        if self.min > len(sequence):
            self.min = len(sequence)


if __name__ == "__main__":
    get_gene_stat(gene_dir, stat, fasta)
