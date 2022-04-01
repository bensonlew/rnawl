# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class RsemDirFile(Directory):
    """
    bam文件夹格式
    """
    def __init__(self):
        super(RsemDirFile, self).__init__()

    def check(self):
        if super(RsemDirFile, self).check():
            return True
