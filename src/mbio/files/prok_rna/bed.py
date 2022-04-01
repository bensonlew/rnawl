# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/26 17:48

import re, os, Bio, argparse, sys, fileinput, urllib2
from biocluster.iofile import *
from  mbio.files.sequence.fasta import *
from mbio.files.gene_structure.gff3 import *
from mbio.files.gene_structure.gtf import *

class BedFile(File):
    
    def __init__(self):
        super(BedFile, self).__init__()
        pass
    
    def check(self):
        super(BedFile,self).check()
        self._check_skechy()
    
    def _check_skechy(self):

        pass
    
    
    
    def to_gtf(self):
        pass
    
    def bed_tbi(self):
        pass
    

