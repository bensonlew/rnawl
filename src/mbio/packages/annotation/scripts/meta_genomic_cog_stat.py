# -*- coding: utf-8 -*-
#!/usr/bin/env python
#author zhouxuan
#last_modify 20170605


import argparse,os,re,collections
import pandas as pd
import numpy

parser = argparse.ArgumentParser(description='to get cog_stat of meta_genomic -- zhouxan')
parser.add_argument("-i", "--cog_table", required=True, help="the result file of the cog_anno")
parser.add_argument("-r", "--reads_profile", required=True, help="the file of the reads_number_profile_table")
parser.add_argument("-o", "--output_dir", required=True, help="the output directory")
args = vars(parser.parse_args())

# 获得cog_anno.xls
func_type = {
	'J': 'INFORMATION STORAGE AND PROCESSING', 'A': 'INFORMATION STORAGE AND PROCESSING',
	'K': 'INFORMATION STORAGE AND PROCESSING', 'L': 'INFORMATION STORAGE AND PROCESSING',
	'B': 'INFORMATION STORAGE AND PROCESSING',
	'D': 'CELLULAR PROCESSES AND SIGNALING', 'Y': 'CELLULAR PROCESSES AND SIGNALING',
	'V': 'CELLULAR PROCESSES AND SIGNALING', 'T': 'CELLULAR PROCESSES AND SIGNALING',
	'M': 'CELLULAR PROCESSES AND SIGNALING', 'N': 'CELLULAR PROCESSES AND SIGNALING',
	'Z': 'CELLULAR PROCESSES AND SIGNALING', 'W': 'CELLULAR PROCESSES AND SIGNALING',
	'U': 'CELLULAR PROCESSES AND SIGNALING', 'O': 'CELLULAR PROCESSES AND SIGNALING',
	'C': 'METABOLISM', 'G': 'METABOLISM', 'E': 'METABOLISM', 'F': 'METABOLISM', 'H': 'METABOLISM',
	'I': 'METABOLISM', 'P': 'METABOLISM', 'Q': 'METABOLISM',
	'R': 'POORLY CHARACTERIZED', 'S': 'POORLY CHARACTERIZED',
}
func_decs = {
	'J': 'Translation, ribosomal structure and biogenesis',
	'A': 'RNA processing and modification', 'K': 'Transcription',
	'L': 'Replication, recombination and repair',
	'B': 'Chromatin structure and dynamics',
	'D': 'Cell cycle control, cell division, chromosome partitioning',
	'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
	'M': 'Cell wall/membrane/envelope biogenesis',
	'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
	'U': 'Intracellular trafficking, secretion, and vesicular transport',
	'O': 'Posttranslational modification, protein turnover, chaperones',
	'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
	'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
	'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
	'P': 'Inorganic ion transport and metabolism',
	'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
	'R': 'General function prediction only', 'S': 'Function unknown'
}

with open(args['cog_table'], 'r') as r, open(args['output_dir'] + "/cog_anno.xls", 'wb') as w:
	w.write("#Query\tCOG/NOG_Group\tFunction\tFunction_des\tFunction_type\n")
	for line in r:
		line = line.strip("\n").split("\t")
		if line[0] != "#Query":
			query_name = line[0]
			OG_Group = line[1]
			function = line[2]
			fun_list = function.split(',')
			fun_d = ('|').join([func_decs[i] for i in fun_list])
			fun_t = ('|').join([func_type[i] for i in fun_list])
			w.write('{}\t{}\t{}\t{}\t{}\n'.format(query_name, OG_Group, function, fun_d, fun_t))

reads_abundance = {}
sample_name = []
with open(args['reads_profile'], 'r') as r:
	for line in r:
		line = line.strip("\n").split("\t")
		if line[0] == "GeneID":
			sample_name = line[1:]
		else:
			reads_abundance[line[0]] = line[1:]

cog_reads = {}
func_reads = {}
cog_des = {}
with open(args['cog_table'], 'r') as r:
	for line in r:
		line = line.strip("\n").split("\t")
		if line[0] != '#Query':
			cog_reads.setdefault(line[1], []).append(line[0])
			func_reads.setdefault(line[2], []).append(line[0])
			if line[1] in cog_des.keys():
				continue
			else:
				if line[4] != '':
					cog_des[line[1]] = line[4]
				else:
					cog_des[line[1]] = line[3]

# 获得cog_list
with open(args['output_dir'] + "/cog_list.xls", 'a') as w:
	w.write('#COG\t' + ('\t').join(sample_name) + '\tDescription\n')
	for key in cog_reads:
		value_list = cog_reads[key]
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
		w.write(key + "\t" + ('\t').join(sum_list) + '\t' + cog_des[key] + '\n')
		# w.write(key + "\t" + ('\t').join(sum_list) + '\tNULL\n')


# 获得dog_func
func_decs = {
	'J': 'Translation, ribosomal structure and biogenesis',
	'A': 'RNA processing and modification', 'K': 'Transcription',
	'L': 'Replication, recombination and repair',
	'B': 'Chromatin structure and dynamics',
	'D': 'Cell cycle control, cell division, chromosome partitioning',
	'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
	'M': 'Cell wall/membrane/envelope biogenesis',
	'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
	'U': 'Intracellular trafficking, secretion, and vesicular transport',
	'O': 'Posttranslational modification, protein turnover, chaperones',
	'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
	'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
	'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
	'P': 'Inorganic ion transport and metabolism',
	'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
	'R': 'General function prediction only', 'S': 'Function unknown'
}
with open(args['output_dir'] + "/cog_function.xls", 'a') as f:
	f.write('#Category\t' + ('\t').join(sample_name) + "\tDescription\n")
	for i in range(65, 91):
		if chr(i) in func_reads.keys():
			reads_list = func_reads[chr(i)]
			if len(reads_list) == 1:
				sum_list = reads_abundance[reads_list[0]]
			else:
				sum_list = []
				for v in reads_list:
					if len(sum_list) == 0:
						sum_list = reads_abundance[v]
					else:
						sum_list = map(int, sum_list)
						reads_abundance[v] = map(int, reads_abundance[v])
						sum_array = numpy.array(sum_list) + numpy.array(reads_abundance[v])
						sum_list = sum_array.tolist()
			sum_list = map(str, sum_list)
			f.write(chr(i) + "\t" + ('\t').join(sum_list) + '\t' + func_decs[chr(i)] + '\n')
		else:
			continue





