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
parser.add_argument('-n', '--num', help='cut num', required=True)
parser.add_argument('-c', '--contain', help='retain length', required=True)
# parser.add_argument('-cut_three', '--cut_three', help='whether to delete the first three bases', required=True)
args = vars(parser.parse_args())

fastq_file = args['fastq']
outfile = args['output']
retain_len = int(args['contain'])
cut_three = "True"
num = int(args['num'])

fs = os.path.basename(fastq_file).split(".")[0] + ".R1.fastq"
if re.search(r".*gz$", fastq_file):
    os.system("gunzip -c {} > {}".format(fastq_file, fs))
else:
    fs = fastq_file
with open(fs, "r") as f, open(outfile, "w") as w:
    i = 1
    for line in f:
        if i % 2 == 0:
            line = line.strip("\n")
            if cut_three == "True":
                w.write(line[num:retain_len+num] + "\n")
            else:
                w.write(line[:retain_len] + "\n")
        else:
            w.write(line)
        i += 1
