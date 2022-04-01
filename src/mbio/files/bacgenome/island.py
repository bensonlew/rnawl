# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""格式文件不为空"""

from biocluster.iofile import File
from biocluster.core.exceptions import OptionError
import os


class IslandFile(File):

    def __init__(self):
        super(IslandFile, self).__init__()

    def check(self):
        if os.path.getsize(self.prop['path']) == 0:
            raise OptionError("文件夹为空", code="41400501")

