# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from biocluster.iofile import File

class CogSummaryFile(File):
    def __init__(self):
        super(CogSummaryFile, self).__init__()

    def check(self):
        if super(CogSummaryFile, self).check():
            return True
