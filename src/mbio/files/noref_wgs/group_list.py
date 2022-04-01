# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class GroupListFile(File):
    """
    wentian@20190110
    """

    def __init__(self):
        super(GroupListFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在%s文件！", variables=(self.path), code="45500307")

    def check(self):
        if super(GroupListFile, self).check():
            self.is_exists()
            with open(self.path, "r") as f:
                lines = f.readlines()
                num = len(lines)
                if num < 2:
                    raise FileError("文件:%s不得少于两行，请检查", variables=(self.path), code="45500308")
                item = lines[0].strip().split("\t")
                if len(item) != 2:
                    raise FileError("文件:%s必须为两列，请检查！", variables=(self.path), code="45500309")
            return True
