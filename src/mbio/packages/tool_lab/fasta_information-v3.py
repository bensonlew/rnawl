#!/usr/bin/env python
# coding: utf8

import re
import argparse

# 读取文件，统计长度
def fasta_infor(fasta_file,out_file):
	dict = {};dict_refer = {}
	seqs_sum = 0;base_sum = 0;GC_sum = 0;N_sum = 0;contig_num = 0
	Ln = 0;base_sum_n = 0
	N50 = 0;N90 = 0
	with open(fasta_file,'r') as read_fasta:
		for line in read_fasta:
			if line[0] == '>':
				key = line.strip('[ >\n]')
				dict[key] = 0
				seqs_sum += 1
			else:
				value = line.strip()
				seqs_len = len(value)
				base_sum += seqs_len
				dict[key] += seqs_len
				GC_sum += (value.count('g') + value.count('G') + value.count('c') + value.count('C'))
				N_sum += (value.count('n') + value.count('N'))
				contig_num += len(re.findall('[nN]+',value))
		contig_num += seqs_sum
		read_fasta.close()
	# 计算N50和N90
	dict_sort = sorted(dict.items(), key = lambda x:x[1], reverse = True)
	base_sum_n50 = 5 * base_sum/10
	base_sum_n90 = 9 * base_sum/10
	for value in dict_sort:
		base_sum_n += value[1]
		Ln += 1
		if N50:
			if base_sum_n >= base_sum_n90:
				N90 = value[1]
				L90 = Ln
				break
		else:
			if base_sum_n >= base_sum_n50:
				N50 = value[1]
				L50 = Ln

	# 输出文件
	average = (base_sum)/(seqs_sum)
	GC =  100 * float(GC_sum) / base_sum
	gap = 100 * float(N_sum) / base_sum
	basic_stat = open(out_file, 'w')
	print >> basic_stat,'Minimum\tLongest\taverage\tGC(%)\tN50\tN90\tgap(%)\tTotal base number\tcontig_num'
	print >> basic_stat,'%d\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d' % (dict_sort[-1][1],dict_sort[0][1],average,GC,N50,N90,gap,base_sum,contig_num)
	basic_stat.close()

def _main():
	parser = argparse.ArgumentParser(description='fasta data')
	parser.add_argument('-i', '--fasta_file', help="fasta_file")
	parser.add_argument('-o', '--out_file', help="out_file")
	args = parser.parse_args()
	fasta_infor(args.fasta_file, args.out_file)

if __name__ == "__main__":
    _main()
