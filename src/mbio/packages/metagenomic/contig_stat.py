# -*- coding: utf-8 -*-
# __author__ = 'guhaidong-20170822'

import os,sys
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description ='file')
#parser.add_argument("-stat_dir", "--stat_dir", required=True)
parser.add_argument("-contig_dir", "--contig_dir", required=True)
#parser.add_argument("-select_stat", "--select_stat", required=False)  #如果没有，不进行kmer筛选
parser.add_argument("-select_kmer", "--select_kmer", required=False)  #是否进行kmer筛选
parser.add_argument("-final_stat", "--final_stat", required=True)
parser.add_argument("-min_contig", "--min_contig", required=False)    #用来统计的最小contig长度
args = vars(parser.parse_args())

contig_dir = args["contig_dir"]
final_stat = args["final_stat"]
final_stat_dir = os.path.dirname(os.path.abspath(final_stat))
if args["select_kmer"] is not None:
    kmer = args["select_kmer"].split(",")
else:
    kmer = None
if args["min_contig"] is None:
    min_contig = 0
else:
    min_contig = int(args["min_contig"])

def get_assemble_stat(contig_dir, final_stat, min_contig, kmer = None):
    all_files = os.listdir(contig_dir)
    name_list = []
    all_stat = dict()
    for files in all_files:
        name_dump = os.path.basename(files).split(".contig")[0]
        one_name = name_dump.split(".kmer")[0]
        if kmer != None:
            one_kmer = name_dump.split(".kmer")[1]
        else:
            if os.path.exists(final_stat_dir + "/" + files):
                sys.stderr.write("相应输出contig文件已存在，跳过此步\n")
            else:
                os.symlink(contig_dir + "/" + files, final_stat_dir + "/" + files)
            try:
                one_kmer = name_dump.split(".kmer")[1]
            except IndexError:
                one_kmer = None
        if one_name not in name_list:
            name_list.append(one_name)
        if one_kmer != None and (kmer != None) and (one_kmer not in kmer):
            raise Exception("错误：路径中存在参数中未设置的kmer！")
        stat_tmp = StatFile(contig_dir + "/" + files, one_kmer, min_contig)
        if one_name not in all_stat:
            all_stat[one_name] = list()
        all_stat[one_name].append(stat_tmp.static)
    with open(final_stat, "w") as fw:
        fw.write("#Sample\tContigs\tContigs bases(bp)\tN50(bp)\tN90(bp)\tMax(bp)\tMin(bp)\n")
        for one in name_list:
            for one_list in all_stat[one]:
                kmer_mark = one_list[-1]
                if kmer == None:
                    stat_line = "\t".join(one_list[0:-1])
                    fw.write(one + "\t" + stat_line + "\n")
                else:
                    [stat_line, the_kmer] = choose_kmer(all_stat[one])
                    stat_line = "\t".join(stat_line)
                    fw.write(one + "\t" + stat_line + "\n")
                    if os.path.exists(final_stat_dir + "/" + one + ".contig.fa"):
                        sys.stderr.write("相应输出contig文件已存在，跳过此步\n")
                    else:
                        os.link(contig_dir + "/" + one + ".kmer" + the_kmer + ".contig.fa", final_stat_dir + "/" + one + ".contig.fa")
                    break
def choose_kmer(stat_list):
    """
    进行kmer筛选
    :param stat_list: 单个样品的contig统计值，二维数组形式存放，每个kmer的统计值为一个数组
    :return:choose_stat,the_kmer 返回最终筛选出的统计值，去掉数组最后一个元素（标记用哪个kmer拼接)
    """
    choose_stat = list()
    the_kmer = ""
    for one in stat_list:
        kmer_tmp = one[-1]
        if len(choose_stat) == 0 or replace(choose_stat, one[0:-1]) == True:
            choose_stat = one[0:-1]
            the_kmer = kmer_tmp
    return choose_stat,the_kmer

def replace(ref_list, replace_list):
    """
    将两个不同kmer的结果进行比较，根据结果决定是否对原最优kmer进行替换
    :param ref_list: 已有的kmer统计结果，目前为最优结果
    :param replace_list: 如果此kmer统计结果比之前的好，则此结果替换掉最优结果
    :return:Boolean  是否替换
    """
    ref_list = map(float,ref_list)
    replace_list = map(float,replace_list)
    if (ref_list[1]/ref_list[0] + ref_list[2]) >= (replace_list[1]/replace_list[0] + replace_list[2]):
        return False
    else:
        return True

class StatFile(object):
    def __init__(self, contig, kmer, min_contig):
        self.contig = contig
        self.n50 = 0
        self.n90 = 0
        self.max = 0
        self.min = 1000000
        self.reads_num = 0
        self.base_num = 0
        self.length_list = []
        self.limit = min_contig
        self.get_stat()
        self.static = [str(self.reads_num), str(self.base_num), str(self.n50), str(self.n90), str(self.max), str(self.min), str(kmer)]

    def get_stat(self):
        """
        对类计算reads数，base数，N50，N90，max_len，min_len
        :return:
        """
        seq = ""
        with open(self.contig) as fr:
            for line in fr:
                line = line.strip()
                if line.startswith(">"):
                    if seq != "":
                        if len(seq) < self.limit:
                            break
                        self.calculate(seq)
                        seq = ""
                else:
                    seq += line
            if len(seq) >= self.limit:
                self.calculate(seq)           #对最后一个序列进行计算
        self.length_list.sort(reverse = True)
        tmp_len = 0
        for element in self.length_list:
            tmp_len += element
            if tmp_len >= self.base_num / 2:
                self.n50 = element
                break
        tmp_len = 0
        for element in self.length_list:
            tmp_len += element
            if tmp_len >= self.base_num * 0.9:
                self.n90 = element
                break

    def calculate(self, sequence):
        """
        计算最大片段，最小片段，总长，总reads数
        :return:
        """
        self.reads_num += 1
        self.base_num += len(sequence)
        self.length_list.append(len(sequence))
        if self.max < len(sequence):
            self.max = len(sequence)
        if self.min > len(sequence):
            self.min = len(sequence)

if __name__ == "__main__":
    get_assemble_stat(contig_dir, final_stat, min_contig, kmer)


