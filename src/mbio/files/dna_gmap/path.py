# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.14

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class PathFile(File):
    """
    检查文件是否存在
    """
    def __init__(self):
        super(PathFile, self).__init__()

    def check(self):
        if super(PathFile, self).check():
            if not os.path.exists(self.path):
                raise FileError("文件:%S不存在，请检查", variables=(self.path), code="44801001")
            else:
                return True


if __name__ == "__main__":
    a = PathFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/01.RawData/fastq.list")
    a.check()
