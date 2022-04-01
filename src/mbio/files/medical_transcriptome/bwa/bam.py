# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class BamFile(File):
    """
    bam文件夹格式
    """
    def __init__(self):
        super(BamFile, self).__init__()

    def check(self):
        if super(BamFile, self).check():
            return True
