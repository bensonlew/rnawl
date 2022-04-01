# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.iofile import File

class GtfFile(File):
    def __init__(self):
        super(GtfFile, self).__init__()

    def check(self):
        super(GtfFile, self).check()
