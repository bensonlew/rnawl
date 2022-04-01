# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""GML格式文件类"""

from biocluster.iofile import File
# import subprocess
# import re
# from biocluster.config import Config
# import os
from biocluster.core.exceptions import FileError


class GmlFile(File):
    """
    Biom文件格式类, 需安装biom工具软件
    """
    def __init__(self):
        """
        """
        super(GmlFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(GmlFile, self).get_info()

    def check(self):
        """
        检测文件是否满足要求
        :return:
        """
        if super(GmlFile, self).check():
            pass
        else:
            raise FileError(u"文件格式错误")
        return True
