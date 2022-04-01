# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class UpsetTableFile(File):
    """
    Binbin Zhao@20200507

    """

    def __init__(self):
        super(UpsetTableFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在{}文件！".format(self.path))

    def check(self):
        if super(UpsetTableFile, self).check():
            self.is_exists()
            # with open(self.path, "r") as f:
            #     lines = f.readlines()
            #     if len(lines[0].strip().split("\t")) != 4:
            #         raise FileError("文件格式有误。数据必须为四列，且以tab分割".format(self.path))
            #     if lines[0].strip().split("\t")[1].isdigit():
            #         raise FileError("文件格式有误。表头不能为空".format(self.path))
            #     for line in lines[1:]:
            #         if not line.strip().split("\t")[1].isdigit():
            #             raise FileError("文件格式有误。第二列必须为数字".format(self.path))
            return True
