# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import File
# from biocluster.core.exceptions import FileError


class CogTableFile(File):
    """
    """
    def __init__(self):
        super(CogTableFile, self).__init__()

    def check(self):
        if super(CogTableFile, self).check():
            return True
