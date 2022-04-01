# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171210

"""fastq序列去掉前三个碱基，保留左边前51bp"""
import argparse
import os
import re


parser = argparse.ArgumentParser(description='Compare the table results to the nt database for species statistics')
parser.add_argument('-i', '--fastq', help='input fastq file', required=True)
parser.add_argument('-o', '--output', help='output stat file', required=True)
# parser.add_argument('-cut_three', '--cut_three', help='whether to delete the first three bases', required=True)
parser.add_argument('-s', '--start', help='start fastq cut', required=True)
parser.add_argument('-l', '--len', help='cut len', required=True)
args = vars(parser.parse_args())

fastq_file = args['fastq']
outfile = args['output']
# cut_three = args['cut_three']
start = int(args['start'])
length = int(args['len'])
end = start + length


fs = os.path.basename(fastq_file).split(".")[0] + ".R1.fastq"
if re.search(r".*gz", fastq_file):
    os.system("gunzip -c {} > {}".format(fastq_file, fs))
else:
    fs = fastq_file
with open(fs, "r") as f, open(outfile, "w") as w:
    i = 1
    for line in f:
        if i % 2 == 0:
            item = line.strip().split()[0]
            line_len = len(item)
            if line_len > end:
                w.write(item[start:end] + "\n")
            elif line_len > start:
                w.write(item[start:] + "\n")
            else:
                w.write(line)
        else:
            w.write(line)
        i += 1
