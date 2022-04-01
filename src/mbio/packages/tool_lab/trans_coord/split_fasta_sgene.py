#!usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import os

parser = argparse.ArgumentParser(description = '\n该脚本用于在基因组特定位置截取序列，需额外输入记录有截取序列信息的列表文件', add_help = False, usage = '\npython3 seq_select.py -i [input.fasta] -o [output.fasta] -l [list]\npython3 seq_select.py --input [input.fasta] --output [output.fasta] --list [list]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = '输入文件，fasta 格式', required = True)
required.add_argument('-o', '--output', metavar = '[output.fasta]', help = '输出文件，fasta 格式', required = True)
args = parser.parse_args()

iput_fasta=args.input
outdir=args.output
n=1
out_prefix='split'
strip_char=".*"
replace=('.', 'X')
with open(iput_fasta) as f:
    seq_num = 1
    chunk_num = 1
    out_name = os.path.join(outdir,"chr" + '_' + str(chunk_num) + '.fa')
    file_objects = {chunk_num: open(out_name, 'w')}
    for line in f:
        if line.startswith("#") or (not line.strip()):
            continue
        if line.startswith('>'):
            chunk_num  = seq_num
            seq_num += 1
            if chunk_num not in file_objects:
                file_objects[chunk_num - 1].close()
                out_name = os.path.join(outdir,"chr" + '_' + str(chunk_num) + '.fa')
                file_objects[chunk_num] = open(out_name, 'w')
        else:
            line = line.strip().strip(strip_char) + '\n'  # 去除两端的*或.
            line = line.replace(replace[0], replace[1])  # 替换'.'为X
        file_objects[chunk_num].write(line)
    else:
        file_objects[chunk_num].close()



