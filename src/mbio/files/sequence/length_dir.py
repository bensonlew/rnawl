# -*- coding: utf-8 -*-
# __author__ = 'sj'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os

class LengthDirFile(File):
    """
    定义gff格式文件
    """

    def __init__(self):
        super(LengthDirFile, self).__init__()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        super(LengthDirFile,self).check()

