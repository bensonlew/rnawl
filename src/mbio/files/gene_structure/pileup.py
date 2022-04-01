# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
# import subprocess


class PileupFile(File):
    """
    定义比对结果sam格式文件
    """

    def __init__(self):
        super(PileupFile, self).__init__()
        # self.samtools_path = "rna/samtools-1.3.1/"

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(PileupFile, self).check():
            return True
        else:
            raise FileError("文件格式错误")