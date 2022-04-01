# -*- coding: utf-8 -*-
# __author__ = 'wentian'
# modified 2019.1.17

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class PopulationsTagFile(File):
    """
    检查文件是否存在
    """
    def __init__(self):
        super(PopulationsTagFile, self).__init__()

    def check(self):
        if super(PopulationsTagFile, self).check():
            if not os.path.exists(self.path):
                raise FileError("文件%s不存在，请检查", variables=(self.path), code="45500905")
            with open(self.path, "r")as fr:
                lines = fr.readlines()
                if len(lines) >= 3:
                    return True
                else:
                    raise FileError("文件%s为空，请检查", variables=(self.path), code="45500906")


if __name__ == "__main__":
    a = PopulationsTagFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/01.RawData/fastq.list")
    a.check()
