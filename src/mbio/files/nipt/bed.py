# -*- coding: utf-8 -*-

# __author__ = 'hongdongxuan'
# time: 2017/5/8 13:29

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class BedFile(File):
    """
    bed文件夹格式
    """
    def __init__(self):
        super(BedFile, self).__init__()

    def check(self):
        if super(BedFile, self).check():
            return True
