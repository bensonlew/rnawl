#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 2021/4/29 16:02
# @Author  : U make me wanna surrender my soul
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", help="bbmap软件路径", type=str, required=True)
parser.add_argument("-r", help="输入生成随机reads的参考序列", type=str, required=True)
parser.add_argument("-o", help="输入生成reads文件的文件名称", type=str, required=True)
parser.add_argument("-l", help="输入生成的reads的长度", type=int, required=True)
parser.add_argument("-paired", help="双端reads or 单端reads", type=str, required=True, default='t')
parser.add_argument("-mininsert", help="输入reads的最小长度", type=int, required=False)
parser.add_argument("-maxinsert", help="输入reads的最大长度", type=int, required=False)
parser.add_argument("-reads", help="输入生成reads的数量", type=int, required=True)
parser.add_argument("-pacbio", help="是否使用Pacbio模型进行reads生成，是 or 否", type=str, required=False, default='f')
parser.add_argument("-split", help="若是双端reads，是否需要分R1，R2,是 or 否", type=str, required=False, default='t')
args = parser.parse_args()
bbmap = args.b
ref = args.r
out = args.o + '.fq.gz'
length = args.l
if args.paired == 't':
        paired = 't'
else:
    paired = 'f'
mininsert = args.mininsert
maxinsert = args.maxinsert
reads = args.reads
if args.pacbio == 'f':
    pacbio = 'f'
else:
    pacbio = 't'
if args.split == 't':
    split_reads = 't'
    cmd_line1 = '{bbmap}/randomreads.sh ref={ref} out={out} length={length} paired={paired} mininsert={mininsert} maxinsert={maxinsert} reads={reads} pacbio={pacbio} overwrite=t'.format(
        bbmap=bbmap, ref=ref, out=out, length=length, paired=paired, mininsert=mininsert, maxinsert=maxinsert, reads=reads,
        pacbio=pacbio)
    print(cmd_line1)
    cmd_line2 = '{bbmap}/reformat.sh in={fastq} out1={R1} out2={R2} addcolon overwrite=t'.format(bbmap=bbmap, fastq=out, R1=args.o + '.R1.fq.gz',
                                                                                     R2=args.o + '.R2.fq.gz')
    subprocess.call(cmd_line1, shell=True)
    subprocess.call(cmd_line2, shell=True)
else:
    split_reads = 'f'
    cmd_line1 = '{bbmap}/randomreads.sh ref={ref} out={out} length={length} paired={paired} mininsert={mininsert} maxinsert={maxinsert} reads={reads} pacbio={pacbio} overwrite=t'.format(
        bbmap=bbmap, ref=ref, out=out, length=length, paired=paired, mininsert=mininsert, maxinsert=maxinsert, reads=reads,
        pacbio=pacbio)
    subprocess.call(cmd_line1, shell=True)
