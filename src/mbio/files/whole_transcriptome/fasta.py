# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import File

class FastaFile(File):
    def __init__(self):
        super(FastaFile, self).__init__()

    def check(self):
        super(FastaFile, self).check()
