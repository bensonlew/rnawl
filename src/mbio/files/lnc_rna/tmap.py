# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class TmapFile(File):
    def __init__(self):
        super(TmapFile, self).__init__()

    def check(self):
        if super(TmapFile, self).check():
            super(TmapFile, self).get_info()
            with open(self.prop["path"], "r") as f:
                for line in f:
                    if line.find("#") == -1:
                        line = line.strip()
                        lst = line.split("\t")
                        if len(lst) != 13:
                            raise FileError("file format error, TMAP should have 13 columns")
                        else:
                            return True
        else:
            raise FileError("file format error")
