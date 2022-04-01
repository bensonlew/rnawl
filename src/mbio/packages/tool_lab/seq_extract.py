#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__ = 'fwy'

#导入模块，初始传递命令、变量等
import argparse

parser = argparse.ArgumentParser(description = '\n该脚本用于在基因组特定位置截取序列，需额外输入记录有截取序列信息的列表文件', add_help = False, usage = '\npython3 seq_select.py -i [input.fasta] -o [output.fasta] -l [list]\npython3 seq_select.py --input [input.fasta] --output [output.fasta] --list [list]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = '输入文件，fasta 格式', required = True)
required.add_argument('-o', '--output', metavar = '[output.fasta]', help = '输出文件，fasta 格式', required = True)
required.add_argument('-l', '--list', metavar = '[list]', help = '记录“新序列名称/序列所在原序列ID/序列起始位置/序列终止位置/正链（+）或负链（-）”的文件，以 tab 作为分隔', required = True)
optional.add_argument('--detail', action = 'store_true', help = '若该参数存在，则在输出 fasta 的每条序列 id 中展示序列在原 fasta 中的位置信息', required = False)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

##读取文件
#读取基因组序列
seq_file = {}
with open(args.input, 'r') as input_fasta:
	for line in input_fasta:
		line = line.strip()
		if line[0] == '>':
			seq_id = line.split()[0]
			seq_file[seq_id] = ''
		else:
			seq_file[seq_id] += line

input_fasta.close()

#读取列表文件
# list_dict = {}
new_name=[]
new_detail=[]
with open(args.list, 'r') as list_file:
	for line in list_file:
		if line.strip():
			line = line.strip().split('\t')
			new_name.append(line[3])
			new_detail.append([line[0], int(line[1]) - 1, int(line[2]), line[4]])

list_file.close()

##截取序列并输出
#定义函数，用于截取反向互补
def rev(seq):
	base_trans = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'a':'t', 'c':'g', 't':'a', 'g':'c'}
	rev_seq = list(reversed(seq))
	rev_seq_list = [base_trans[k] for k in rev_seq]
	rev_seq = ''.join(rev_seq_list)
	return(rev_seq)

#截取序列并输出
output_fasta = open(args.output, 'w')
for n,name in enumerate(new_name):
	if args.detail:
		output_fasta.write('>' + name+ ' '+ str(new_detail[n][0])+' '+ str(new_detail[n][1] + 1) +' '+str(new_detail[n][2])+' '+str(new_detail[n][3]) + '\n')
	else:
		output_fasta.write('>' + name)
	
	seq = seq_file['>' + new_detail[n][0]][new_detail[n][1]:new_detail[n][2]]
	if new_detail[n][3] == '+':
		output_fasta.write(seq+"\n")
	elif new_detail[n][3] == '-':
		seq = rev(seq)
		output_fasta.write(seq+"\n")

output_fasta.close()
