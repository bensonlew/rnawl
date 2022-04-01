# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re,shutil
import argparse
from Bio import SeqIO

def get_gff3(samples, faa_dir, gff_dir,genome_dir, out):
    """
    根据输入的样本名称，蛋白文件夹、gff文件夹、全基因组文件夹，输出文件夹
    """
    if os.path.exists(out):
        shutil.rmtree(out)
    os.mkdir(out)
    #print samples
    all_sample_list = samples.strip().split("+")
    all_sample_list.sort()
    for sample in all_sample_list:
        #print sample
        sample_faa_name = sample + '_CDS.faa'
        sample_gff_name = sample + '_CDS.gff'
        faa_path = os.path.join(faa_dir, sample_faa_name)
        gff_path = os.path.join(gff_dir, sample_gff_name)
        out_gff = os.path.join(out, sample_gff_name)
        merge_table(sample, faa_path, gff_path, out_gff)

def merge_table(sample, faa, gff, out_gff):
    """
    针对单个样本合并生成gff文件
    :param sample:
    :param faa:
    :param gff:
    :param out_gff:
    :return:
    """
    with open(gff, 'r') as f, open(out_gff, 'w') as w:
        for line in f:

            w.write(line)
        w.write("##FASTA\n")
        for seq_record in SeqIO.parse(faa, 'fasta'):  # 打开第一次主要是获取到样本名称
            seq_id = seq_record.id
            seq_name = sample + "|" + seq_id
            seq = seq_record.seq
            w.write(">{}\n{}\n".format(seq_name, seq))

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[sample]',required=True, help='input sample name')
    parse.add_argument('-faa', metavar='[faa_dir]', required=True, help='input faa dir')
    parse.add_argument('-gff', metavar='[gff_dir]', required=True, help='input gff dir')
    parse.add_argument('-genome', metavar='[genome_dir]', required=True, help='input genome dir')
    parse.add_argument('-o', metavar='[output_dir]', required=True, help='input output dir')
    args = parse.parse_args()
    samples = args.i
    faa_dir = args.faa
    gff_dir = args.gff
    genome_dir = args.genome
    out = args.o
    get_gff3(samples,faa_dir, gff_dir, genome_dir, out)

if __name__ == '__main__':
    main()
