# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import File

class CommonFile(File):
    def __init__(self):
        super(CommonFile, self).__init__()

    def check(self):
        super(CommonFile, self).check()
