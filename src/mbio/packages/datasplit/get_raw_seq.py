# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171214

"""获取多样性的原始fastq"""
import os
import re
import argparse


parser = argparse.ArgumentParser(description='获取多样性的原始fastq')
parser.add_argument('-r1', '--fastq1', help='*.raw.valid.1.fq', required=True)
parser.add_argument('-r2', '--fastq2', help='*.raw.valid.2.fq', required=True)
parser.add_argument('-s', '--seq2sam', help='*.raw.seq2sam.stat', required=True)
parser.add_argument('-o', '--outdir', help='outdir', required=True)

args = vars(parser.parse_args())
fastq1 = args['fastq1']
fastq2 = args['fastq2']
seq2sam = args['seq2sam']
outdir = args['outdir']

seq_sample = {}
sample_list = []
with open(seq2sam, "rb") as f:
    lines = f.readlines()
    for line in lines[1:]:
        item = line.strip().split("\t")
        seq_sample[item[0]] = item[1]
        if item[1] not in sample_list:
            sample_list.append(item[1])
            # break
        # print item[0]
# print seq_sample

r1_files = {}
r2_files = {}
for s in sample_list:
    r1_files[s] = open(os.path.join(outdir, s + ".R1.raw.fastq"), "w")
    r2_files[s] = open(os.path.join(outdir, s + ".R2.raw.fastq"), "w")

with open(fastq1, "rb") as fastq:
    for line in fastq:
        seq_id = line.split(" ")[0]
        # print seq_id
        if seq_id in seq_sample:
            s = seq_sample[seq_id]
            r1_files[s].write("{}{}{}{}".format(line, next(fastq), next(fastq), next(fastq)))
        # break
with open(fastq2, "rb") as fastq:
    for line in fastq:
        seq_id = line.split(" ")[0]
        # print seq_id
        if seq_id in seq_sample:
            s = seq_sample[seq_id]
            r2_files[s].write("{}{}{}{}".format(line, next(fastq), next(fastq), next(fastq)))
        # break

for s in sample_list:
    r1_files[s].close()
    r2_files[s].close()
