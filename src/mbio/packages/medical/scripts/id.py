# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
# last_modify:2016.12.5

import sys
import re
# s = WQ2131-FC

# id_modified.py /mnt/ilustre/users/sanger-dev/workspace/20161205/Single_fastq2tab_module/Fastq2tab/output/bam2tab WQ2131-M
# path = sys.argv[1]
file = sys.argv[1]

with open(file,'r+') as f:
	for line in f:
		line = line.strip()
		line = line.split('\t')
		line[0] = line[0][:-1]
		dst = '\t'.join(line)
		s = re.match(".*/output/(.*)\.tab", file)
		x = s.group(1)
		tab = open(x + "_rename.tab","a+")
		tab.write(dst+'\n')
		tab.close()
