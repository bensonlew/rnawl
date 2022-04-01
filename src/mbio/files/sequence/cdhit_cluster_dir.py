# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import re
import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class CdhitClusterDirFile(Directory):
    """
    定义CdihitClustr文件夹
    """

    def __init__(self):
        super(CdhitClusterDirFile, self).__init__()

    def get_info(self):
        """
        获取文件夹属性
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            self.set_property("files_number", self.get_number())
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", code="44000401")

    def get_number(self):
        filelist = os.listdir(self.prop['path'])
        count = 0
        count_dirs = 0
        for file_ in filelist:
            if os.path.isfile(os.path.join(self.prop['path'], file_)):
                count += 1
            if os.path.isdir(os.path.join(self.prop['path'], file_)):
                count_dirs += 1
        return count_dirs

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(CdhitClusterDirFile, self).check():
            self.get_info()
            return True
