# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""文件夹类"""

from biocluster.iofile import File
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import re
import os


class PathDirFile(Directory):
    def __init__(self):
        super(PathDirFile, self).__init__()

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(PathDirFile, self).check():
            if not os.path.isdir(self.prop["path"]):
                raise FileError("文件格式错误")
        else:
            raise FileError("文件格式错误")
