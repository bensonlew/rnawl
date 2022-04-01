#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kefei.huang
# 20171121

import pysam
import re

samfile = pysam.AlignmentFile("WQ17114164-F1.mem.sort.hit.bam","rb")
outfile = pysam.AlignmentFile("WQ17114164-F1.mem.sort.hit.bam.filter.bam","wb",template=samfile)
for line in samfile:
	if line.alen < 10:
		continue
	if line.mapping_quality < 30:
		continue
	if re.search(r"S",line.cigarstring) != None:
		continue
	outfile.write(line)

samfile.close()
outfile.close()
	