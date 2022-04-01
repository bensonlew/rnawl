# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import collections
from Bio import SeqIO

class ArfToWig(object):
    def __init__(self):
        self.arf = ""
        self.wig = ""
        self.step_width = 5000
        self.mapping = collections.OrderedDict()
        self.fa = ""
        self.sample = ""
        self.indexs = []
        self.chr_list = []

    def set_arf(self, arf):
        '''
        设置arf路径
        '''
        self.arf = arf

    def set_fa(self, fa):
        '''
        设置fa路径
        '''
        self.fa = fa

    def init_mapping(self):
        '''
        初始化列表
        '''
        index = 0
        for seq in SeqIO.parse(self.fa, "fasta"):
            chrom = seq.id
            length = len(seq.seq)
            start = 1
            while True:
                end = start + self.step_width - 1
                if end > length:
                    end = length
                    pos = (start + end)/2
                    self.mapping.update({chrom + "_" + str(start): {
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "pos": pos,
                        "num_pos": 0,
                        "num_neg": 0
                    }})
                    break
                else:
                    pos = (start + end)/2
                    self.mapping.update({chrom + "_" + str(start): {
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "pos": pos,
                        "num_pos": 0,
                        "num_neg": 0
                    }})
                    start = end + 1

    def parser_arf(self):
        '''
        统计arf文件
        '''
        with open(self.arf, 'rb') as f:
            print self.arf
            for line in f:
                cols = line.strip().split("\t")
                sam = cols[0].split("_")[0]
                if self.sample != "" and sam != self.sample:
                    continue
                num = int(cols[0].split("_")[2].lstrip("x"))
                start = int(cols[7])
                end = int(cols[8])
                chrom = cols[5]
                length  = int(cols[1])
                start_pre = start/self.step_width * self.step_width + 1
                index_pre = chrom + '_' + str(start_pre)
                index_next = chrom + '_' + str(start_pre + self.step_width)
                if cols[10] == "+":
                    num_strand = "num_pos"
                else:
                    num_strand = "num_neg"
                if end >= start_pre + self.step_width:
                    length_n = end - (start_pre + self.step_width) + 1
                    length_p = length - length_n
                    self.mapping[index_pre][num_strand] = self.mapping[index_pre][num_strand] + length_p/length * num
                    self.mapping[index_next][num_strand] = self.mapping[index_next][num_strand] + length_n/length * num
                else:
                    self.mapping[index_pre][num_strand] = self.mapping[index_pre][num_strand] + num

    def write_file(self, stat_file):
        '''
        输出统计相关结果
        '''
        with open(stat_file, 'wb') as f:
            for index,values in self.mapping.items():
                f.write("\t".join([values['chr'], str(values['start']), str(values['end']), str(values['pos']), str(values['num_pos']), str(values['num_neg'])]) + '\n')

    def write_window(self, stat_file):
        '''
        输出统计相关结果
        '''
        with open(stat_file + '.pos.window', 'wb') as f_pos, open(stat_file + '.neg.window', 'wb') as f_neg:
            pos_list = [value_dict['num_pos'] for value_dict in self.mapping.values()]
            neg_list = [value_dict['num_neg'] for value_dict in self.mapping.values()]
            num_max = max(pos_list + neg_list)
            chr_set = set(self.chr_list)
            for index,values in self.mapping.items():
                if values['chr'] in chr_set:
                    f_pos.write("\t".join([values['chr'], str(values['start']), str(values['end'] - 1), str(float(values['num_pos'])/num_max), "fill_color=col6"])  + '\n')
                    f_neg.write("\t".join([values['chr'], str(values['start']), str(values['end'] - 1), str(float(values['num_neg'])/num_max), "fill_color=col3"])  + '\n')


if __name__ == "__main__":
    test_dir = '/mnt/ilustre/users/sanger-dev/workspace/20181112/Single_MapperAndStat14-35-13/MapperAndStat/output/'
    A = ArfToWig()
    A.step_width = 500000
    A.set_fa(test_dir + 'Caenorhabditis_elegans.WBcel235.dna_rm.toplevel.fa')
    A.set_arf(test_dir + 'reads_vs_genome.arf')
    A.init_mapping()
    # print A.mapping
    A.parser_arf()
    A.chr_list=['I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA']
    A.sample = 'S02'
    A.write_file('stat.xls')
    A.write_window('S02_circos')
