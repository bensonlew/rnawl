# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

""".gz检查"""
from biocluster.iofile import File
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError
import re


class FileGzFile(File):
    """
    txt类
    """

    def __init__(self):
        super(FileGzFile, self).__init__()

    def check(self):
        """
    检测文件是否以.gz结尾
    :return:
        """
        if not re.search(r'.gz',self.prop['path']):
            raise FileError('文件格式不是以.gz结尾的文件!')
        else:
            pass
        return True