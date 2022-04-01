# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20210730

"""从fastq里随机抽取num条组成fasta"""
import os
import re
import random
import argparse


parser = argparse.ArgumentParser(description='从fastq里随机抽取num条组成fasta')
parser.add_argument('-fq1', '--fq1', help='fq1', required=True)
parser.add_argument('-fq2', '--fq2', help='fq2', required=False)
parser.add_argument('-num', '--num', help='num', required=True)
parser.add_argument('-out', '--out', help='out', required=True)

args = parser.parse_args()
fq1 = args.fq1
if fq1.endswith(".gz"):
    new_fq = os.path.basename(fq1).split(".gz")[0]
    os.system("gunzip -c {} > {}".format(fq1, new_fq))
    fq1 = new_fq
if args.fq2:
    fq2 = args.fq2
    if fq2.endswith(".gz"):
        new_fq = os.path.basename(fq2).split(".gz")[0]
        os.system("gunzip -c {} > {}".format(fq2, new_fq))
        fq2 = new_fq
num = int(args.num)
out = args.out

def get_rand_list(path, rand_num):
    count = 0
    thefile = open(path, 'rb')
    while True:
        buffer = thefile.read(8192*1024)
        if not buffer:
            break
        count += buffer.count('\n')
    thefile.close()
    print "文件行数:%s" % count
    seq_num = count / 4
    num_dict = {}
    if seq_num <= rand_num:
        for i in range(seq_num):
            num_dict[i] = 1
    else:
        while len(num_dict.keys()) < rand_num and len(num_dict.keys()) < seq_num:
            rand = random.randint(0, seq_num)
            if rand not in num_dict:
                num_dict[rand] = 1
    print "从0-%s里随机抽取%s个数" % (seq_num, len(num_dict.keys()))
    return num_dict

def get_new_fasta(path, num_dict, name, out_fa):
    count, i = 0, 0
    num = len(num_dict.keys())
    with open(path) as fastq:
        for line in fastq:
            i = i + 1
            seq = next(fastq)
            plus = next(fastq)
            qual = next(fastq)
            if i in num_dict:
                out_fa.write(">"+name+str(count)+"\n" + seq)
                count = count + 1
            if count == num:
                break

rand_num = num
if args.fq2:
    rand_num = num / 2
num_dict1 = get_rand_list(fq1, rand_num)
name = os.path.basename(fq1).split("-")[-1] + "1_"
out_fa = open(out, "wb")
get_new_fasta(fq1, num_dict1, name, out_fa)
if args.fq2:
    num_dict2 = get_rand_list(fq2, rand_num)
    name = os.path.basename(fq2).split("-")[-1] + "2_"
    get_new_fasta(fq2, num_dict2, name, out_fa)
out_fa.close()
