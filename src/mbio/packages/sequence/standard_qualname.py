# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# __version__ = 'v1.0'
# __last_modified__ = '20160112'
"""
检查fasta序列名和qual序列名的第一空位符前的名称是否一致，如果一致，将qual的名称转换成fasta的名称，否则报错
"""
import os
# import argparse
# import sys


# def get_argu():
#     """
#     命令行模式下运行此脚本的参数获取方法
#     """
#     __parse_details = argparse.ArgumentParser(prog='standard_qualname.py',
#                                     usage='关于此脚本的说明',
#                                     description='检查fasta序列名和qual序列名的第一空位符前的名称是否一致，\
#                                     如果一致，将qual的名称转换成fasta的名称，否则报错',
#                                     epilog='请确保参数准确',
#                                     version='v1.0',
#                                     parents='')
#     __parse_details.add_argument('-f', '--fasta',
#                        required=True, help='输入fasta文件')
#     __parse_details.add_argument('-q', '--qual',
#                        required=True, help='输入fasta对应qual文件')
#     __parse_details.add_argument('-o', '--output',
#                        required=True, help='输出新的qual文件路径')
#     args = __parse_details.parse_args()
#
#     output = args.output
#     fasta = args.fasta
#     qual = args.qual
#     return fasta, qual, output


def standard_qualname(fasta, qual, newqual):
    def set_error():
        fastafile.close()
        qualfile.close()
        new.close()
        os.remove(newqual)

    try:
        fastafile = open(fasta)
    except IOError:
        return 1, 'fasta文件(%s)无法打开' % fasta
    try:
        qualfile = open(qual)
    except IOError:
        return 1, 'qual质量文件(%s)无法打开' % qual
    try:
        new = open(newqual, 'w')
    except IOError:
        return 1, '无法创建新的newqual文件(%s)' % newqual
    while True:
        fastaname = ''
        qualname = ''
        fastaend = False
        qualend = False
        while True:
            fastaline = fastafile.readline()
            if fastaline == '':
                fastaend = True
                break
            elif fastaline[0] == '>':
                fastaname = fastaline
                break
            else:
                pass

        while True:
            qualline = qualfile.readline()
            if qualline == '':
                qualend = True
                break
            elif qualline[0] == '>':
                qualname = qualline
                break
            else:
                new.write(qualline)
        if not qualend and not fastaend:
            pass
        elif qualend and fastaend:
            return 0, '正常'
        elif not fastaend and qualend:
            set_error()
            return 1, 'fasta 与 qual 序列量不一样'
        elif not qualend and fastaend:
            set_error()
            return 1, 'fasta 与 qual 序列量不一样'
        else:
            return 1, '未知错误'
        if fastaname.split()[0] == qualname.split()[0]:
            new.write(fastaname)
        else:
            return 1, 'fasta名称(%s)与qual名称(%s)不一致' % (fastaname.split()[0], qualname.split()[0])

# if __name__ == '__main__':
#     fasta, qual, newqual = get_argu()
#     returncode = standard_qualname(fasta, qual, newqual)
#     if returncode[0] == 0:
#         sys.exit(0)
#     else:
#         sys.exit(1)
