#!/usr/bin/env python
#author zhouxuan
#last_modify 20170602

import argparse,os,re,collections
import pandas as pd
import numpy

parser = argparse.ArgumentParser(description='to get nr_tax_profile of on meta_genomic -- zhouxan')
parser.add_argument("-i", "--tax_file", required=True, help="the result file of the nr_anno")
parser.add_argument("-r", "--reads_profile", required=True, help="the file of the reads_number_profile_table")
parser.add_argument("-o", "--output_dir", required=True, help="the output directory")
args = vars(parser.parse_args())

dict_1 = {}
with open(args['tax_file'], 'r') as r:
	for line in r:
		line = line.strip("\n").split("\t")
		dict_1.setdefault(line[1], []).append(line[0])

dict_2 = {}
sample_name= []
with open(args['reads_profile'], 'r') as r:
	for line in r:
		line = line.strip("\n").split("\t")
		if line[0] == "GeneID":
			sample_name = line[1:]
		else:
			dict_2[line[0]] = line[1:]

with open(args['output_dir'] + "/tax_profile.xls", 'a') as w:
	first_line = "#Taxonomy\t" + ('\t').join(sample_name) + "\n"
	w.write(first_line)
	for key in dict_1:
		value_list = dict_1[key]
		if len(value_list) == 1:
			sum_list = dict_2[value_list[0]]
		else:
			sum_list = []
			for v in value_list:
				if len(sum_list) == 0:
					sum_list = dict_2[v]
				else:
					sum_list = map(int, sum_list)
					dict_2[v] = map(int, dict_2[v])
					sum_array = numpy.array(sum_list) + numpy.array(dict_2[v])
					sum_list = sum_array.tolist()
		sum_list = map(str, sum_list)
		w.write(key + "\t" + ('\t').join(sum_list) + '\n')



