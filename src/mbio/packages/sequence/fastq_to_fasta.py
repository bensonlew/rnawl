# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# __version__ = 'v1.0'
# __last_modified__ = '20160112'
"""
将fastq文件转化为fasta文件以及其qual文件
"""


import os
# import argparse
# import sys


# def get_argu():
#     """
#     命令行模式下运行此脚本的参数获取方法
#     """
#     __parse_details = argparse.ArgumentParser(prog='fastq_to_fasta.py',
#                                     usage='关于此脚本的说明',
#                                     description='将fastq文件转换成fasta文件，并可以生成qual质量文件，phred默认为33',
#                                     epilog='请确保参数准确',
#                                     version='v1.0',
#                                     parents='')
#     __parse_details.add_argument('-f', '--fastq',
#                        required=True, help='输入fastq文件')
#     __parse_details.add_argument('-a', '--fasta',
#                        required=True, help='输出的fasta文件')
#     __parse_details.add_argument('-o', '--qual',
#                        required=True, help='输出qual质量文件')
#     __parse_details.add_argument('-p', '--phred',
#                        required=True, help='质量的格式')
#     args = __parse_details.parse_args()
#
#     fastq = args.fastq
#     fasta = args.fasta
#     qual = args.qual
#     phred = args.phred
#     return fastq, fasta, qual, phred


def convertfastq(fastq, outfasta, outqual='', phred=33):
    def set_error():
        pass
    try:
        fastqfile = open(fastq)
    except IOError:
        return 1, '无法打开fastq文件：%s' % fastq
    try:
        fastafile = open(outfasta, 'w')
    except IOError:
        return 1, '无法创建fasta文件：%s' % outfasta
    if outqual != '':
        try:
            qualfile = open(outqual, 'w')
        except IOError:
            return 1, '无法创建qual文件：%s' % outqual
    flag = 0
    line1 = ''
    line2 = ''
    line3 = ''
    line4 = ''
    for line in fastqfile:
        if flag == 0:
            line1 = line
        elif flag == 1:
            line2 = line
        elif flag == 2:
            line3 = line
        elif flag == 3:
            print 'flag:', flag
            line4 = line
            flag = -1
            if line1[0] == '@' and line3[0] == '+':
                fastafile.write('>' + line1[1:])
                fastafile.write(line2)
                if outqual:
                    qualfile.write('>' + line1[1:])
                    qualfile.write(ord_seq(line4, phred))
            else:
                return 1, '错误的fastq格式：%s%s%s%s' % (line1, line2, line3, line4)
        else:
            return 1, '未知错误'
        flag = flag + 1


def ord_seq(seq, phred):
    new_seq = ''
    for i in seq.strip():
        new_seq += (str(ord(i) - phred) + ' ')
    return new_seq.strip() + '\n'

# if __name__ == '__main__':
#     fastq, fasta, qual, phred = get_argu()
#     return_code = convertfastq(fastq, fasta, qual, phred)
#     if return_code[0] == None:
#         sys.exit(0)
#     elif return_code[0] == 1:
#         sys.exit(1)
#         print return_code[1]
#     else:
#         print '未知错误'
#         sys.exit(2)
