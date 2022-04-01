# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
##程序是对sam处理，
##1. 删除没有序列信息的
##2. 删除比对到了softcliping的
##3. 删除mapping quality < 30的

import codecs
import re
from __future__ import print_function

SAMOUT = codecs.open("WQ17114164-F1.sam_filter.sam","w","utf-8")
with codecs.open("WQ17114164-F1.sam","r","utf-8") as SAM:
	for line in SAM:
		line = line.strip()
		if line.startswith("@"):
			print (line,file=SAMOUT)
		else:
			linetemp = line.split("\t")
			if (re.search(r"S",linetemp[5])) != None:
				print (line,file=ERROROUT)
				continue
			if (linetemp[4] < 30):
				print (line,file=ERROROUT)
				continue
			if (len(linetemp[9]) < 2):
				print (line,file=ERROROUT)
				continue
			print (line,file=SAMOUT)
SAMOUT.close()