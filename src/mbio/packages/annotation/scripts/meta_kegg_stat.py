# -*- coding: utf-8 -*-
#!/usr/bin/env python
#author zhouxuan
#last_modify 20170608


import argparse, os, re, collections
import pandas as pd
import numpy

parser = argparse.ArgumentParser(description='to get kegg_stat of meta_genomic -- zhouxan')
parser.add_argument("-k", "--kegg_table", required=True, help="the result file of the kegg_anno")
parser.add_argument("-e", "--enzyme_table", required=True, help="the result file of the enzy_table")
parser.add_argument("-p", "--path_table", required=True, help="the result file of the path_table")
parser.add_argument("-m", "--module_table", required=True, help="the result file of the module_table")
parser.add_argument("-r", "--reads_profile", required=True, help="the file of the reads_number_profile_table")
parser.add_argument("-o", "--output_dir", required=True, help="the output directory")
args = vars(parser.parse_args())

# 获得gene.profile和KO.profile文件
gene = {}
KO = {}
gene_KO = {}
with open(args['kegg_table'], 'r') as k:
	for line in k:
		line = line.strip("\n").split("\t")
		if line[0] != "#Query":
			gene.setdefault(line[1], []).append(line[0])
			KO.setdefault(line[2], []).append(line[0])
			gene_KO[line[0]] = line[2]  # 后面pathway的时候用
		else:
			continue
reads_abundance = {}
sample_name = []
with open(args['reads_profile'], 'r') as r:
	for line in r:
		line = line.strip("\n").split("\t")
		if line[0] == "GeneID":
			sample_name = line[1:]
		else:
			reads_abundance[line[0]] = line[1:]
with open(args['output_dir'] + "/gene_profile", 'wb') as g:
	g.write('gene\t' + ('\t').join(sample_name) + '\n')
	for key in gene:
		value_list = gene[key]
		if len(value_list) == 1:
			sum_list = reads_abundance[value_list[0]]
		else:
			sum_list = []
			for v in value_list:
				if len(sum_list) == 0:
					sum_list = reads_abundance[v]
				else:
					sum_list = map(int, sum_list)
					reads_abundance[v] = map(int, reads_abundance[v])
					sum_array = numpy.array(sum_list) + numpy.array(reads_abundance[v])
					sum_list = sum_array.tolist()
		sum_list = map(str, sum_list)
		g.write(key + "\t" + ('\t').join(sum_list) + '\n')
with open(args['output_dir'] + "/KO_profile", 'wb') as K:
	K.write('KO\t' + ('\t').join(sample_name) + '\n')
	for key in KO:
		value_list = KO[key]
		if len(value_list) == 1:
			sum_list = reads_abundance[value_list[0]]
		else:
			sum_list = []
			for v in value_list:
				if len(sum_list) == 0:
					sum_list = reads_abundance[v]
				else:
					sum_list = map(int, sum_list)
					reads_abundance[v] = map(int, reads_abundance[v])
					sum_array = numpy.array(sum_list) + numpy.array(reads_abundance[v])
					sum_list = sum_array.tolist()
		sum_list = map(str, sum_list)
		K.write(key + "\t" + ('\t').join(sum_list) + '\n')

# Module.profile文件的生成
module_des = {}
module = {}
with open(args["module_table"], 'r') as m:
	for line in m:
		line = line.strip("\n").split("\t")
		module_list = line[1].split(",")
		module_dlist = line[2].split(";")
		t = len(module_list)
		for i in range(0, t):
			if module_list[i] not in module_des.keys():
				module_des[module_list[i]] = module_dlist[i]
			else:
				continue
		for M in module_list:
			module.setdefault(M, []).append(line[0])
with open(args['output_dir'] + "/module_profile.xls", 'wb') as m:
	m.write('Module\t' + ('\t').join(sample_name) + '\tDefinition\n')
	for key in module:
		value_list = module[key]
		if len(value_list) == 1:
			sum_list = reads_abundance[value_list[0]]
		else:
			sum_list = []
			for v in value_list:
				if len(sum_list) == 0:
					sum_list = reads_abundance[v]
				else:
					sum_list = map(int, sum_list)
					reads_abundance[v] = map(int, reads_abundance[v])
					sum_array = numpy.array(sum_list) + numpy.array(reads_abundance[v])
					sum_list = sum_array.tolist()
		sum_list = map(str, sum_list)
		m.write(key + "\t" + ('\t').join(sum_list) + "\t" + module_des[key] + '\n')

# enzyme.profile生成
enzyme_des = {}
enzyme = {}
with open(args["enzyme_table"], 'r') as e:
	for line in e:
		line = line.strip("\n").split("\t")
		enzyme_list = line[1].split(",")
		enzyme_dlist = line[2].split(";")
		t = len(enzyme_list)
		for i in range(0, t):
			if enzyme_list[i] not in enzyme_des.keys():
				enzyme_des[enzyme_list[i]] = enzyme_dlist[i]
			else:
				continue
		for M in enzyme_list:
			enzyme.setdefault(M, []).append(line[0])
with open(args['output_dir'] + "/enzyme_profile.xls", 'wb') as e:
	e.write('Enzyme\t' + ('\t').join(sample_name) + '\tDefinition\n')
	for key in enzyme:
		value_list = enzyme[key]
		if len(value_list) == 1:
			sum_list = reads_abundance[value_list[0]]
		else:
			sum_list = []
			for v in value_list:
				if len(sum_list) == 0:
					sum_list = reads_abundance[v]
				else:
					sum_list = map(int, sum_list)
					reads_abundance[v] = map(int, reads_abundance[v])
					sum_array = numpy.array(sum_list) + numpy.array(reads_abundance[v])
					sum_list = sum_array.tolist()
		sum_list = map(str, sum_list)
		e.write(key + "\t" + ('\t').join(sum_list) + "\t" + enzyme_des[key] + '\n')

# pathway.profile生成
path_way_des = {}
path_way = {}
with open(args["path_table"], 'r') as p:
	for line in p:
		line = line.strip("\n").split("\t")
		path_list = line[1].split(",")
		path_dlist = line[2].split(";")
		t = len(path_list)
		for i in range(0, t):
			if path_list[i] not in path_way_des.keys():
				path_way_des[path_list[i]] = path_dlist[i]
			else:
				continue
		for M in path_list:
			path_way.setdefault(M, []).append(line[0])
with open(args['output_dir'] + "/pathway_profile.xls", 'wb') as p:
	p.write('Pathway\t' + ('\t').join(sample_name) + '\tDefinition\tPathwaymap\n')
	for key in path_way:
		value_list = path_way[key]
		if len(value_list) == 1:
			sum_list = reads_abundance[value_list[0]]
		else:
			sum_list = []
			for v in value_list:
				if len(sum_list) == 0:
					sum_list = reads_abundance[v]
				else:
					sum_list = map(int, sum_list)
					reads_abundance[v] = map(int, reads_abundance[v])
					sum_array = numpy.array(sum_list) + numpy.array(reads_abundance[v])
					sum_list = sum_array.tolist()
		sum_list = map(str, sum_list)
		ko_list = []
		for i in value_list:
			if gene_KO[i] in ko_list:
				pass
			else:
				ko_list.append(gene_KO[i])
		path_map = ("+").join(ko_list)
		path_m = 'http://www.genome.jp/kegg-bin/show_pathway?' + key + '+' + path_map
		p.write(key + "\t" + ('\t').join(sum_list) + "\t" + path_way_des[key] + "\t" + path_m + '\n')