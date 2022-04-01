# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

"""metagenome原始数据统计运行时需要的文件夹"""

from biocluster.iofile import Directory
from biocluster.core.exceptions import OptionError
import os


class RawdataDirFile(Directory):
    """
    特定形式的fastq.gz文件夹
    """
    def __init__(self):
        super(RawdataDirFile, self).__init__()
        
    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        allfiledir = os.listdir(self.prop['path'])
        if len(allfiledir) == 0:
            raise OptionError("文件夹为空")
        else:
            return True
