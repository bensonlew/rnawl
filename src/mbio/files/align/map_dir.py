# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class MapDirFile(Directory):
    """
    bwt_index文件夹格式
    """
    def __init__(self):
        super(MapDirFile, self).__init__()

    def check(self):
        if super(MapDirFile, self).check():
            return True