# -*- coding: utf-8 -*-

# __author__ = 'moli.zhou'
# time: 2017/5/8 13:29

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class XlsxFile(File):
    """
    bed文件夹格式
    """
    def __init__(self):
        super(XlsxFile, self).__init__()

    def check(self):
        if super(XlsxFile, self).check():
            return True
