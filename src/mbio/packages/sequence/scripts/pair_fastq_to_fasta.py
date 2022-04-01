#!/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# __version__ = 'v1.0'
# __last_modified__ = '20160112'
"""
将成对fastq文件转换成fasta文件
"""

import argparse


def get_argu():
    par = argparse.ArgumentParser()
    par.add_argument("-f", metavar="[fastqfile1]", required=True, help="输入成对fastq文件1")
    par.add_argument("-r", metavar="[fastqfile2]", required=True, help="输入成对fastq文件2")
    par.add_argument("-o", metavar="[outfile]", required=True, help="输出fasta文件")
    par.add_argument("-n", metavar="[id_name1]", type=str, default="none", help="可自定义fastq1转换成fasta文件的id名字")
    par.add_argument("-m", metavar="[id_name2]", type=str, default="none", help="可自定义fastq1转换成fasta文件的id名字")
    args = par.parse_args()
    return args


def pair_fastq_to_fasta(fq1, fq2, fasta, id1, id2):
    """
    将fastq文件转换成fasta文件
    :param fq1:成对fastq文件1
    :param fq2:成对fastq文件2
    :return:
    """
    n = 1
    with open(fq1, 'r') as r1:
        with open(fq2, 'r') as r2:
            with open(fasta, 'w') as w:
                file1 = r1.readlines()
                file2 = r2.readlines()
                length = len(file1)
                for i in range(1, length+1):
                    if (i-1) % 4 == 0:
                        if id1 == "none" and id2 == "none":
                            w.write('%s' % file1[i-1].replace('@', '>'))
                            w.write(file1[i])
                            w.write('%s' % file2[i-1].replace('@', '>'))
                            w.write(file2[i])
                        elif id1 != "none" and id2 != "none":
                            w.write('%s%s\n' % (id1, n))
                            w.write(file1[i])
                            w.write('%s%s\n' % (id2, n))
                            w.write(file2[i])
                            n += 1
                        elif id1 != "none" and id2 == "none":
                            w.write('%s%s\n' % (id1, n))
                            w.write(file1[i])
                            w.write('%s' % file2[i-1].replace('@', '>'))
                            w.write(file2[i])
                            n += 1
                        else:
                            w.write('%s' % file1[i-1].replace('@', '>'))
                            w.write(file1[i])
                            w.write('%s%s\n' % (id2, n))
                            w.write(file2[i])
                            n += 1


def main():
    opts = get_argu()
    pair_fastq_to_fasta(opts.f, opts.r, opts.o, opts.n, opts.m)
if __name__ == '__main__':
    main()
