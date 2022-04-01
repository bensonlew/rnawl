# -*- coding: utf-8 -*-
#!/usr/bin/env python
#author zhouxuan
#last_modify 20170612


import argparse, os, re, collections
import pandas as pd
import numpy

parser = argparse.ArgumentParser(description='to get nr_level of meta_genomic -- zhouxan')
parser.add_argument("-i", "--nr_tax_table", required=True, help="the result file of the nr_stat")
parser.add_argument("-l", "--level", required=False, help="the level you want", default="1,2,3,4,5,6,7,8")
parser.add_argument("-o", "--output_dir", required=True, help="the output directory")
args = vars(parser.parse_args())

level_dict = {'1': 'd',
              '2': 'k',
              '3': 'p',
              '4': 'c',
              '5': 'o',
              '6': 'f',
              '7': 'g',
              '8': 's'}
# 获取需要获得结果表的分类级别
level = args['level'].split(',')
for i in level:
	level_number = {}
	sample_name = []
	with open(args['nr_tax_table'], 'r') as n:
		for line in n:
			if line.startswith("#Taxonomy"):
				sample_name = line.strip('\n').split('\t')[1:]
			else:
				line = line.strip('\n').split('\t')
				level_name = (';').join(line[0].split(';')[0:int(i)])
				level_number.setdefault(level_name, []).append(line[1:])
	with open(args['output_dir'] + "/tax_" + level_dict[i] + ".xls", 'wb') as t:
		t.write('#Taxonomy\t' + ('\t').join(sample_name) + '\n')
		for key in level_number:
			value_list = level_number[key]
			if len(value_list) == 1:
				sum_list = level_number[key][0]
			else:
				sum_list = []
				for v in value_list:
					if len(sum_list) == 0:
						sum_list = v
					else:
						sum_list = map(int, sum_list)
						v = map(int, v)
						sum_array = numpy.array(sum_list) + numpy.array(v)
						sum_list = sum_array.tolist()
			sum_list = map(str, sum_list)
			t.write(key + "\t" + ('\t').join(sum_list) + '\n')

#  python metagen_nr_taxlevel.py -i /mnt/ilustre/users/sanger-dev/workspace/20170612/Single_nr_stat_6251/MetagenNr
# /output/tax_profile.xls -l 1,2,3,4,5,6,7 -o /mnt/ilustre/users/sanger-dev/sg-users/zhouxuan/metagenome/nr_taxlevel
