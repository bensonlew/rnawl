# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""bam类"""

from biocluster.iofile import Directory
from biocluster.core.exceptions import OptionError


class SampleDataDirFile(Directory):
    """
    bam文件夹
    """
    def __init__(self):
        super(SampleDataDirFile, self).__init__()
        
    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(SampleDataDirFile, self).check():
            return True
        else:
            raise FileError("文件格式错误")
