# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/16 13:29

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class AsFile(File):
    """
    as文件夹格式
    """
    def __init__(self):
        super(AsFile, self).__init__()

    def check(self):
        if super(AsFile, self).check():
            return True