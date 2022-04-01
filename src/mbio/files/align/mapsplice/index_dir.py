# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/25 10:00

import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class IndexDirFile(Directory):
    """
    mapsplice bowtie index文件夹格式
    """
    def __init__(self):
        super(IndexDirFile, self).__init__()

    def check(self):
        if super(IndexDirFile, self).check():
            return True