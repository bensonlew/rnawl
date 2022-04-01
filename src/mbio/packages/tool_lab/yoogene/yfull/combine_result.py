#!/bin/env python
# coding=utf-8

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputtable", help="antosome结果表")
parser.add_argument('-n',"--name", help="sample_name")
args = vars(parser.parse_args())

with open(args['inputtable'], 'r') as f, open(args['name']+'_antosome_snp_result.txt', 'w') as w:
	lines = f.readlines()
	w_sep = '\t'
	line_dict = {}
	result_dict = {}
	rs_list = []
	last_list = []
	title_list = lines[0].strip().split('\t')
	w.write(w_sep.join(title_list))
	w.write('\n')
	for line in lines[1:]:
		line_list = line.strip().split('\t')
		if line_list[0] in line_dict.keys():
			rs_list = line_dict[line_list[0]]
		else:
			pass
		if line_list[2] not in rs_list:
			if line_list[0] not in result_dict.keys():
				result_dict[line_list[0]] = {}
			result_dict[line_list[0]][line_list[2]] = line_list
			line_dict[line_list[0]] = line_list[2]
		else:
			if line_list[0] in result_dict.keys():
				if line_list[2] in result_dict[line_list[0]].keys():
					last_list = result_dict[line_list[0]][line_list[2]]
			l_list0 = ''
			l_list1 = ''
			n_str0 = ''
			n_str1 = ''
			if len(last_list) != 0:
				if '/' in list(last_list[3]):
					l_list0 = last_list[3].split('/')[0]
					l_list1 = last_list[3].split('/')[1]
				else:
					for n in range(len(list(last_list[3]))):
						if n < len(list(last_list[3]))/2:
							l_list0 += list(last_list[3])[n]
						else:
							l_list1 += list(last_list[3])[n]
				for m in range(len(list(line_list[3]))):
					if m < len(list(line_list[3]))/2:
						n_str0 += list(line_list[3])[m]
					else:
						n_str1 += list(line_list[3])[m]
				new_str = l_list0 + n_str0 + '/' + l_list1 + n_str1
				last_list[3] = new_str
				last_pos = last_list[1]
				last_list[1] = last_pos + ";" + line_list[1]
				last_DP = last_list[4]
				last_list[4] = last_DP + ";" + line_list[4]
				result_dict[line_list[0]][line_list[2]] = last_list
	for v in result_dict.values():
		w_list = v.values()
		for t in w_list:
			w.write(w_sep.join(t))
			w.write('\n')
