# coding=utf-8

from Bio import SeqIO
import fileinput
import re
import os
import subprocess
import urllib2
from collections import defaultdict
import regex
# import pandas as pd
# from collections import Counter
# import matplotlib.pyplot as plt
import sys

def step_count(fasta_file, fasta_to_txt, group_num, step, stat_out, min_len=0):
    """
    步长统计
    :param fasta_file: 输入的fa文件
    :param fasta_to_txt:输出的统计数据txt
    :param group_num:按照步长统计几组
    :param step:步长
    :param stat_out:统计的数据汇总信息txt（fasta_to_txt文件的汇总）
    :return:
    """
    with open(fasta_to_txt, "w") as f:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            ID = seq_record.description.strip().split(" ")[0]
            new_trans_list = ID + "\t" + str(len(seq_record)) + "\n"
            f.write(new_trans_list)
    with open(fasta_to_txt, "r") as r, open(stat_out, "a") as w:
        sample_name = os.path.basename(fasta_file).split('.fa')[0]
        w.write("Length" + "\tNumber" + "\n")
        trans_list = []
        amount_group = []
        element_set = set("")
        for line in r:
            line = line.strip().split("\t")
            number = line[1]
            trans_list.append(number)
        for f in trans_list:
            for i in range(group_num):
                if (int(f) >= (i * step)) and (int(f) < ((i+1) * step)):
                    amount_group.append(i)
                else:
                    pass
                    # amount_group.append(group_num+1)
                element_set.add(i)
        amount_group.sort()
        top_sum = 0
        all_sum = sum([amount_group.count(i) for i in element_set])

        for i in element_set:
            num_statistics = amount_group.count(i)
            # pct = format(float(num_statistics)/all_sum, '.2%')
            if str(i) == '0':
                area_line = str(min_len) + "~" + str((i + 1) * step) + "\t" + str(num_statistics) +  "\n"
                w.write(area_line)
                top_sum += int(num_statistics)
            elif i < (group_num-1):
                area_line = str(i * step + 1) + "~" + str((i+1) * step) + "\t" + str(num_statistics)  + "\n"
                w.write(area_line)
                top_sum += int(num_statistics)
            else:
                area_line = ">" + str(i * step) + "\t" + str(len(trans_list)-int(top_sum))+ "\n"
                w.write(area_line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', type=str, metavar="fasta_file", required=True,
                        help="please input fasta file ")
    parser.add_argument('-s', type=int, metavar="step", default=200, help="step for split length ", required=True)
    parser.add_argument('-n', type=int, metavar="group_num", default=10, help="group number to split", required=True)
    parser.add_argument('-o', type=str, metavar="output_dir",default=None, help="default is local dir. Output directory name", required=True)
    parser.add_argument('-m', type=str, metavar="min_length", default=None, help="min length for stat")
    #
    args = parser.parse_args()
    file_path,out,group_num,step  = args.f ,args.o , args.n ,args.s
    step_count(file_path,os.path.join(out,"seq_length_deatail"),group_num,step,os.path.join(out,"seq_statistic"),0)
