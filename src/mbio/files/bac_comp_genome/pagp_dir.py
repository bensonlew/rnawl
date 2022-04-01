# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2019.05.10

from biocluster.iofile import Directory
import os,re
from biocluster.core.exceptions import FileError

class PagpDirFile(Directory):
    """
    PAGP输入文件夹的检查
    """
    def __init__(self):
        super(PagpDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(PagpDirFile, self).get_info()
        file_list = os.listdir(self.prop['path'])
        for file in  file_list:
            if not file.endswith((".pep", ".nuc", ".function")):
                raise FileError('%s文件后缀格式不正确，请检查！', variables=(file))

    def check(self):
        if super(PagpDirFile, self).check():
            self.get_info()
            return True

if __name__ == "__main__":
    dir = PagpDirFile()