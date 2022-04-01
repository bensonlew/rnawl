# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class KeggListFile(File):
    """
    定义kegg.list格式
    """
    def __init__(self):
        super(KeggListFile, self).__init__()

    def check(self):
        if super(KeggListFile, self).check():
            return True
