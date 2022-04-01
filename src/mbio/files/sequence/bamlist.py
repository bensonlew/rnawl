# -*- coding: utf-8 -*-

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import subprocess
import os

class BamlistFile(File):
    """
    定义Bamlist格式文件
    """

    def __init__(self):
        super(BamlistFile, self).__init__()

    def get_info(self):
        super(BamlistFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        super(BamlistFile,self).check()
