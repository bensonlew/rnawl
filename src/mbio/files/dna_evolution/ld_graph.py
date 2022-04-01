# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.08.22

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class LdGraphFile(File):
    """
    定义windowed.pi格式
    必须是五列，不能为空文件
    """
    def __init__(self):
        super(LdGraphFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在{}文件！".format(self.path))

    def check(self):
        if super(LdGraphFile, self).check():
            self.is_exists()
            return True