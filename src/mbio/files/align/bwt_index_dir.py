# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class BwtIndexDirFile(Directory):
    """
    bwt_index文件夹格式
    """
    def __init__(self):
        super(BwtIndexDirFile, self).__init__()

    def check(self):
        if super(BwtIndexDirFile, self).check():
            return True