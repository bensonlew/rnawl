# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

import os
import re
import argparse


def seqstat(fasta):
    cmd = '/mnt/ilustre/users/sanger/app/bioinfo/seqs/biosquid_1.9g+cvs20050121/bin/seqstat -a %s  >> seqstat_outfile' % fasta
    os.system(cmd)
    record = open('seqstat_outfile', 'r')
    out_file = open('fasta_len.xls', 'w')
    out_file_sec = open('fasta_stat.xls', 'w')
    lines = record.readlines()
    for l in lines[-7:]:
        out_file_sec.write('%s' % l)
    out_file.write('\n')
    # print lines[5]
    for i in re.findall('[\S]+', lines[5]):
        # print i
        out_file.write('%s\t' % i)
    out_file.write('\n')
    for l in lines[7:-8]:
        # print l
        for i in re.findall('[\S]+', l)[1:]:
            out_file.write('%s\t' % i)
        out_file.write('\n')


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', '--infile', help='input fasta file')
    args = parse.parse_args()
    _i = args.infile
    seqstat(_i)

if __name__ == '__main__':
    main()
