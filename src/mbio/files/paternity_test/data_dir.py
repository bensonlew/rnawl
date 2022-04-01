# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

"""亲子鉴定数据拆分输入文件夹运行时需要的data文件夹"""

from biocluster.iofile import Directory
from biocluster.core.exceptions import OptionError
import os


class DataDirFile(Directory):
    """
    特定形式的fastq.gz文件夹
    """
    def __init__(self):
        super(DataDirFile, self).__init__()
        
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
