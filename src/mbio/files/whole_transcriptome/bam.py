# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import File

class BamFile(File):
    def __init__(self):
        super(BamFile, self).__init__()

    def check(self):
        super(BamFile, self).check()
