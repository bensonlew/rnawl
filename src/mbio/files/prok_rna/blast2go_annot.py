# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import File
# from biocluster.core.exceptions import FileError


class Blast2goAnnotFile(File):
    """
    """
    def __init__(self):
        super(Blast2goAnnotFile, self).__init__()

    def check(self):
        if super(Blast2goAnnotFile, self).check():
            return True
