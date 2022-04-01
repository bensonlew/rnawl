# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

from biocluster.iofile import File

class BamFile(File):
    def __init__(self):
        super(BamFile, self).__init__()

    def check(self):
        if super(BamFile, self).check():
            return True
