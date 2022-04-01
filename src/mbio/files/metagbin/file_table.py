# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""amphora注释"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError


class FileTableFile(File):
    """
    txt类
    """

    def __init__(self):
        super(FileTableFile, self).__init__()

    def check(self):
        """
    检测文件是否为空，数据是否存在
    :return:
        """
        with open(self.prop['path'], 'r') as r:
            lines = r.readlines()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(FileTableFile, self).get_info()