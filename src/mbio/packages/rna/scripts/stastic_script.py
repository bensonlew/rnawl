## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
import sys

file = sys.argv[1]
row = int(sys.argv[2]) - 1
print row
s = dict()

# file = 'TF_result_plant.txt'
# row =1
with open(file,'r+') as f:
	for line in f:
		line = line.strip()
		line = line.split('\t')
		temp = line[row]
		if s.has_key(temp):
			s[temp] = int(s[temp]) + 1
		else:
			s[temp] = int(1)


for key in s:
	fs = open('stastic_result.txt', 'a+')
	fs.write(key+'\t'+str(s[key])+'\n')
	fs.close()





