# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
import os
import re
import gzip
import subprocess
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fasta import FastaFile

class GeneInFastaFile(FastaFile):
	def __init__(self):
		super(GeneInFastaFile, self).__init__()

	def check_gene(self):
		with open(self.prop['path'], 'r') as r:
			line = r.next()
			if re.search(r'>.*', line):
				line = line.strip()
				g = re.search('gene:([a-zA-Z0-9\.]*)\s', line)
				gene = g.group(1)
				if not gene:
					raise FileError('文件中不包含差异基因对应信息')
		return True