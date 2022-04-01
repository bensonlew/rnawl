# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import File

class FastqFile(File):
    def __init__(self):
        super(FastqFile, self).__init__()

    def check(self):
        super(FastqFile, self).check()
