# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import Directory


class DiffStatDirFile(Directory):
    """
    差异基因文件夹格式
    """
    def __init__(self):
        super(DiffStatDirFile, self).__init__()

    def check(self):
        if super(DiffStatDirFile, self).check():
            return True
