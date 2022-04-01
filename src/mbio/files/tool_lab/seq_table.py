# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class SeqTableFile(File):
    """
    wuqin @20200908
    """

    def __init__(self):
        super(SeqTableFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在{}文件！".format(self.path))

    def check(self):
        if super(SeqTableFile, self).check():
            self.is_exists()
            with open(self.path, "r") as f:
                lines = f.readlines()
                if len(lines) != 1:
                    raise FileError("{}文件格式有误。序列写在文件第一行，不能换行".format(self.path))
                for line in lines:
                    seq = str(line.strip().split("\t")[0])
                    if " " in seq or "," in seq or "." in seq:
                        raise FileError("{}文件格式有误。输入序列中不能有空格或','或'.'或Tab键空格".format(self.path))
            return True
