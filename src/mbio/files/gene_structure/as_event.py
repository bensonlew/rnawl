# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/26 17:48

import re, os, Bio, argparse, sys, fileinput, urllib2
from biocluster.iofile import *
from  mbio.files.sequence.fasta import *
from mbio.files.gene_structure.gff3 import *
from mbio.files.gene_structure.gtf import *


class AsEventFile(File):
    def __init__(self):
        super(AsEventFile, self).__init__()
        pass

    def check(self):
        """
            检查文件是否满足要求，不满足触发FileError异常
        """
        super(AsEventFile, self).check()
        if super(AsEventFile, self).check():
            return True
        else:
            raise FileError("as_event文件格式错误！")

