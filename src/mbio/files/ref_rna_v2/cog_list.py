# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import File
# from biocluster.core.exceptions import FileError


class CogListFile(File):
    """
    """
    def __init__(self):
        super(CogListFile, self).__init__()

    def check(self):
        if super(CogListFile, self).check():
            return True
