# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import Directory

class CommonDirFile(Directory):
    def __init__(self):
        super(CommonDirFile, self).__init__()

    def check(self):
        super(CommonDirFile, self).check()
