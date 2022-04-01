## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
import sys
import os

file = sys.argv[1]
path = sys.argv[2]
# file = 'Ca_HPV.RDDpred.results_report.txt'

with open(file,'r+') as f:
	for line in f:
		line = line.strip()
		line = line.split('\t')
		ref_base = line[2]
		var_base = line[3]
		ratio = line[5]
		if ratio == 'Prediction_Likelihood_Ratio':
			continue
		elif float(ratio) <= 0.1:
			temp = '0-0.1'
		elif 0.1 < float(ratio) <= 0.2:
			temp = '0.1-0.2'
		elif 0.2 < float(ratio) <= 0.3:
			temp = '0.2-0.3'
		elif 0.3 < float(ratio) <= 0.4:
			temp = '0.3-0.4'
		elif 0.4 < float(ratio) <= 0.5:
			temp = '0.4-0.5'
		elif 0.5 < float(ratio) <= 0.6:
			temp = '0.5-0.6'
		elif 0.6 < float(ratio) <= 0.7:
			temp = '0.6-0.7'
		elif 0.7 < float(ratio) <= 0.8:
			temp = '0.7-0.8'
		elif 0.8 < float(ratio) <= 0.9:
			temp = '0.8-0.9'
		elif 0.9 < float(ratio) <= 1:
			temp = '0.9-1'

		fs = open(os.path.join(path, 'edit_reassin.txt'), "a+")
		# fs = open('edit_reassin.txt', "a+")
		fs.write(ref_base+'2'+var_base+'\t'+temp+'\n')
		fs.close





