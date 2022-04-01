# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class VcfDirFile(Directory):
    """
    定义vcf文件
    """
    def __init__(self):
        super(VcfDirFile, self).__init__()

    def check(self):
        """
        检查文件是否满足要求，不满足触发FileError异常
        """
        if super(VcfDirFile, self).check():
            return True
        else:
            raise FileError("文件格式错误！")
    """

     def get_info(self):

        super(VcfFile, self).get_info()
        self.get_vcf_info()



    def get_vcf_info(self):
        with open(self.prop["path"],"r") as v:
            lines = v.readlines()

    """