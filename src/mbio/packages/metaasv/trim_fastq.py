# -*- coding: utf-8 -*-

import sys
import pandas as pd
import argparse
import logging
import os
import shutil
from Bio import SeqIO


def find_min_length(file):
    """
    找出fastq序列文件的最小长度
    :param file:
    :return:
    """
    min_length = []
    for seq_record in SeqIO.parse(file, "fastq"):
        seq = str(seq_record.seq)
        length = len(seq)
        if length not in min_length:
            min_length.append(length)
    min_seq_number = min(min_length)
    return min_seq_number


def trim_data(input_table, output_table, trim_length):
    """
    对输入序列进行标准化，挑选固定长度的序列
    :param input_table: 输入fastq序列文件
    :param output_table: 输出fastq序列文件
    :param trim_length: trim_length
    :return:
    """
    with open(input_table, "r") as f, open(output_table, "w") as w:
        for line in f:
            if line[0] == "@":
                first_line = line
                seq = next(f)
                if len(seq) > int(trim_length):
                    w.write(first_line)
                    w.write(seq[0:int(trim_length)]+"\n")
                    w.write(next(f))
                    quality = next(f)
                    w.write(quality[0:int(trim_length)]+"\n")


def _main(input, output, min_number=True):
    """
    首先求出每个样本的最小序列长度
    然后根据最小序列长度过滤数据
    :param input: 输入文件夹
    :param output: 输出文件夹
    :param min_number :最小样本数
    :return:
    """
    all_list = os.listdir(input)
    if not min_number:
        min_seq = []
        for file in all_list:
            file_path = os.path.join(input, file)
            seq_number = find_min_length(file_path)
            if seq_number not in min_seq:
                min_seq.append(seq_number)
        total_min_number = min(min_seq)
    else:
        total_min_number = min_number
    if total_min_number < 90:
        total_min_number = 90
    print("real_min_number：{}".format(total_min_number))
    for f in all_list:
        old_file = os.path.join(input, f)
        new_file = os.path.join(output,f)
        trim_data(old_file, new_file, total_min_number)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[fastq_dir]',required=True,help='fastq dir')
    parser.add_argument('-o', metavar='[output_dir]', required=True, help='output_dir')
    parser.add_argument('-min', metavar='[min_sample_number]', help='min sample number')
    args = parser.parse_args()
    input_dir = args.i
    out_dir = args.o
    if args.min:
        min_number = args.min
        _main(input_dir, out_dir, min_number)
    else:
        _main(input_dir, out_dir)

