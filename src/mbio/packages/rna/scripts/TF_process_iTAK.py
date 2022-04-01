# -*- coding: utf-8 -*-
# author: Moli
# last update: 20161118

import re
import sys
import os
import subprocess

# python TF_process_iTAK.py /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/ /mnt/ilustre/users/sanger-dev/app/bioinfo/rna/iTAK-1.6b/ Arabidopsis_thaliana.TAIR10.pep.all.fa /mnt/ilustre/users/sanger-dev/workspace/20161116/Single_express_sample_v8_change_sample_name/Express/output/diff/featurecounts_count.txt.A_vs_B.edgeR.DE_results_name

cmd = '{}perl {}iTAK.pl {}' \
	.format(sys.argv[1], sys.argv[2], sys.argv[3])

try:
	subprocess.check_output(cmd, shell=True)
# res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError:
	raise Exception("运行出错！")

id = []
family = []
TR = []
result_path = sys.argv[3] + '_output'
tf_clasification = result_path + "/" + 'tf_classification.txt'
with open(tf_clasification, 'r') as f:
	for line in f:
		line = line.strip()
		line = line.split('\t')
		id.append(line[0])
		family.append(line[1])
		TR.append(line[2])


i = -1
pep_id = []
gene_id = []
seq = []
with open(sys.argv[3], 'r+') as f:
	for line in f:
		line = line.strip()
		if '>' in line:
			x = re.match('>([a-zA-Z0-9\.]*)\s', line)
			pep_id.append(x.group(1))
			y = re.search('gene:([a-zA-Z0-9\.]*)\s', line)
			gene_id.append(y.group(1))
			i += 1
			seq.append('')
		else:
			seq[i] += line + '\n'

diff = []
with open(sys.argv[4], 'r') as d:
	for line in d:
		line = line.strip()
		diff.append(line)

for k in range(len(id)):
	if family[k] == "Other" or TR[k] != "TF":
		pass
	else:
		for m in range(len(pep_id)):
			if id[k] == pep_id[m]:
				if gene_id[m] in diff:
					f_out = open("TF_result.txt", "a+")
					f_out.write(pep_id[k] + '\t' + family[k] + '\t' + gene_id[k] + '\n')
					f_out.close()
				else:
					f_out = open("TF_result.txt", "a+")
					f_out.write(pep_id[k] + '\t' + family[k] + '\t' + 'no' + '\n')
					f_out.close()

