# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import re
import os
import sys
import math
import pandas as pd


class BedIntersect(object):
    def __init__(self):
        '''
        self.bed_protein 存储蛋白位置信息
        self.bed_lnc 存储lncRNA位置信息
        self.length_max 判断是否为靶基因最大长度
        self.length_bidirection 判断是否为bidirectional_promoter_lncRNA 的最大长度
        '''

        self.bed_protein = dict()
        self.bed_lnc = dict()
        self.length_max = 10000
        self.length_max_up = 10000
        self.length_max_down = 10000
        self.length_bidirection = 1000
        self.lnctype = dict()
        self.lnc_types = ['intergenic', 'unstrand_overlap', "unstrand_intron_overlap", "unstrand_exon_overlap", 'bidirection', "antisense", "antisense_intron_overlap", "antisense_exon_overlap", "sense_intron_overlap", "sense_exon_overlap"]
        self.windows_mrna = dict()

    def get_lnc_pos(self, lnc_bed):
        '''
        由bed文件获取lnc信息, bed应该为bed12格式
        '''
        with open(lnc_bed, 'rb') as f:
            for line in f:
                cols = line.strip("\n").split("\t")
                lnc_id = cols[3]
                lnc_chr = cols[0]
                lnc_pos1 = int(cols[1])
                lnc_pos2 = int(cols[2])
                lnc_strand = cols[5]
                exons_size = cols[10].strip(",").split(",")
                exons_start = cols[11].strip(",").split(",")
                exons = list()
                for i,s in enumerate(exons_size):
                    exons.append([lnc_pos1 + int(exons_start[i]) + 1,  lnc_pos1 + int(exons_start[i]) + int(s)])

                self.bed_lnc[lnc_id] = [lnc_chr, lnc_pos1, lnc_pos2, lnc_strand, exons]


    def get_protein_pos(self, protein_bed):
        '''
        由bed文件获取protein信息, bed应该为bed12格式
        '''
        with open(protein_bed, 'rb') as f:
            for line in f:
                cols = line.strip("\n").split("\t")
                protein_id = cols[3]
                protein_chr = cols[0]
                protein_pos1 = int(cols[1])
                protein_pos2 = int(cols[2])
                protein_strand = cols[5]
                exons_size = cols[10].strip(",").split(",")
                exons_start = cols[11].strip(",").split(",")
                exons = list()
                for i,s in enumerate(exons_size):
                    exons.append([protein_pos1 + int(exons_start[i]) + 1,  protein_pos1 + int(exons_start[i]) + int(s)])
                self.bed_protein[protein_id] = [protein_chr, protein_pos1, protein_pos2, protein_strand, exons]
                '''
                建立搜索滑窗，目的是减少搜索范围
                '''
                window_start = protein_pos1/100000
                window_end = protein_pos2/100000
                for win in range(window_start, window_end + 1):
                    win_key = protein_chr + "_" + str(win)
                    if win_key in self.windows_mrna:
                        self.windows_mrna[win_key].add(protein_id)
                    else:
                        self.windows_mrna[win_key] = set([protein_id])

        # print self.windows_mrna

    def check_exon_overlap(self, exons1, exons2):
        '''
        获取外显子是否交叉
        '''
        for range1 in exons1:
            for range2 in exons2:
                if range1[1] >= range2[0] and range1[0] >= range2[1]:
                    return True
        return False

    def init_lnctype(self, lnctype):
        '''
        初始化lncRNA类型
        '''
        for lnc_id in self.bed_lnc.keys():
            self.lnctype[lnc_id] = lnctype

    def change_lnctype(self, lnc_id, new_type):
        '''
        修改lncRNA类型, 根据类型顺序用后面的替换前面的类型
        lnc_types = ['intergenic', 'bidirection', "antisense_intron_overlap", "antisense_exon_overlap", "sense_intron_overlap", "sense_exon_overlap"]
        '''
        if self.lnc_types.index(new_type) > self.lnc_types.index(self.lnctype[lnc_id]):
            self.lnctype[lnc_id] = new_type

    def lnc_type(self, out_put):
        '''
        获取靶基因类型
        '''
        self.init_lnctype("intergenic")

        for lnc_id, lnc_ele in self.bed_lnc.items():
            search_start = lnc_ele[1] - self.length_bidirection
            search_end = lnc_ele[2] + self.length_bidirection
            search_chr = lnc_ele[0]
            window_start = int(search_start)/100000
            window_end = int(search_end)/100000
            proteins = set()
            for win in range(window_start, window_end + 1):
                win_key = search_chr + "_" + str(win)
                # print lnc_id, win_key
                if win_key in self.windows_mrna:
                    proteins = proteins | self.windows_mrna[win_key]
            # print proteins
            # for protein_id, protein_ele in self.bed_protein.items():
            for protein_id in proteins:
                protein_ele = self.bed_protein[protein_id]
                if protein_ele[0] == lnc_ele[0]:
                    inter_type = "intergenic"
                    if lnc_ele[3] == "+" and (protein_ele[1] - lnc_ele[2]) > 0:
                        continue
                    elif lnc_ele[3] == "+":
                        if lnc_ele[1] - protein_ele[2] > self.length_bidirection:
                            continue
                        else:
                            if protein_ele[3] == "-":
                                inter_type = "bidirection"
                    elif lnc_ele[3] == "-" and lnc_ele[1] - protein_ele[2] > 0:
                        continue
                    elif lnc_ele[3] == "-":
                        if protein_ele[1] - lnc_ele[2] > self.length_bidirection:
                            continue
                        else:
                            if protein_ele[3] == "+":
                                inter_type = "bidirection"
                    if protein_ele[1] - lnc_ele[2] <= 0 and lnc_ele[1] - protein_ele[2] <= 0:
                        if protein_ele[3] == "." or  lnc_ele[3] == ".":
                            inter_type = "unstrand_overlap"
                        elif protein_ele[3] == lnc_ele[3]:
                            inter_type = "sense"
                            if self.check_exon_overlap(protein_ele[4], lnc_ele[4]):
                                inter_type += "_exon_overlap"
                            else:
                                inter_type += "_intron_overlap"
                        else:
                            inter_type = "antisense"

                    # print inter_type, lnc_id, protein_id
                    if inter_type != "intergenic":
                        self.change_lnctype(lnc_id, inter_type)

        with open(out_put, 'w') as f:
            for lnc_id,lnc_type in self.lnctype.items():
                f.write("{}\t{}\n".format(lnc_id, lnc_type))

        types = self.lnctype.values()
        with open(out_put + ".stat.xls", 'w') as f:
            for lnc_type in self.lnc_types:
                f.write("{}\t{}\n".format(lnc_type, types.count(lnc_type)))

    def intersection(self, out_put):
        '''
        获取邻近区域靶基因信息
        '''
        with open(out_put, 'w') as f:
            f.write("lncRNA\tmRNA\tlncRNA_pos\tlncRNA_strand\tmRNA_pos\tmRNA_strand\ttype\tdistance\n")
            for protein_id, protein_ele in self.bed_protein.items():
                for lnc_id, lnc_ele in self.bed_lnc.items():
                    if protein_ele[0] == lnc_ele[0]:
                        if lnc_ele[3] == "+" and (protein_ele[1] - lnc_ele[2] > self.length_max_down or lnc_ele[1] - protein_ele[2] > self.length_max_up):
                            pass
                        elif lnc_ele[3] == "-" and (protein_ele[1] - lnc_ele[2] > self.length_max_up or lnc_ele[1] - protein_ele[2] > self.length_max_down):
                            pass
                        else:
                            inter_type = "intergenic"
                            if protein_ele[1] - lnc_ele[2] > 0:
                                distance = protein_ele[1] - lnc_ele[2]
                                if lnc_ele[3] == '+':
                                    inter_type += "_downstream"
                                if lnc_ele[3] == '-':
                                    inter_type += "_upstream"
                            elif lnc_ele[1] - protein_ele[2] > 0:
                                distance = lnc_ele[1] - protein_ele[2]
                                if lnc_ele[3] == '-':
                                    inter_type += "_downstream"
                                if lnc_ele[3] == '+':
                                    inter_type += "_upstream"
                            else:
                                distance = 0

                            if protein_ele[1] - lnc_ele[2] > 0 and protein_ele[1] - lnc_ele[2] < self.length_bidirection:
                                if protein_ele[3] == "+" and lnc_ele[3] == "-":
                                    inter_type = "bidirection"
                            if lnc_ele[1] - protein_ele[2] > 0 and lnc_ele[1] - protein_ele[2] < self.length_bidirection:
                                if protein_ele[3] == "-" and lnc_ele[3] == "+":
                                    inter_type = "bidirection"

                            if protein_ele[1] - lnc_ele[2] <= 0 and lnc_ele[1] - protein_ele[2] <= 0:
                                if protein_ele[3] == "." or  lnc_ele[3] == ".":
                                    inter_type = "unstrand"
                                elif protein_ele[3] == lnc_ele[3]:
                                    inter_type = "sense"
                                else:
                                    inter_type = "antisense"

                                if self.check_exon_overlap(protein_ele[4], lnc_ele[4]):
                                    inter_type += "_exon_overlap"
                                else:
                                    inter_type += "_intron_overlap"

                            f.write("\t".join([lnc_id, protein_id, protein_ele[0],
                                               str(lnc_ele[1] + 1) + "-" + str(lnc_ele[2]),
                                               lnc_ele[3],
                                               str(protein_ele[1] + 1) + "-" + str(protein_ele[2]),
                                               protein_ele[3],
                                               inter_type,
                                               str(distance)
                            ]) + "\n")

if __name__ == "__main__":
    BedInter = BedIntersect()
    if (len(sys.argv) == 4):
        BedInter.get_lnc_pos(sys.argv[1])
        BedInter.get_protein_pos(sys.argv[2])
        BedInter.lnc_type(sys.argv[3])
    elif len(sys.argv) == 6 :
        BedInter.length_max_up = int(sys.argv[4]) * 1000
        BedInter.length_max_down = int(sys.argv[5]) * 1000
        BedInter.get_lnc_pos(sys.argv[1])
        BedInter.get_protein_pos(sys.argv[2])
        BedInter.intersection(sys.argv[3])
