# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.06.29

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class MarkerListFile(File):
    """
    检查文件是否存在
    """
    def __init__(self):
        super(MarkerListFile, self).__init__()

    def check(self):
        if super(MarkerListFile, self).check():
            if not os.path.exists(self.path):
                raise FileError("文件:%s不存在，请检查", variables=(self.path), code="44800501")
            self.check_column()

    def check_column(self):
        """
        检查文件列数是否为2列。
        """
        with open(self.path, "r")as fc:
            lines = fc.readlines()
            for line in lines:
                marker_list = line.strip().split("\t")
                list_len = len(marker_list)
                if list_len != 2:
                    raise FileError("文件格式不对，不是两列！", code="44800502")

if __name__ == "__main__":
    a = MarkerListFile()
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/gmap/05.map-cycle1-18-1-4/ref.marker.list")
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/gmap/05.map-cycle1-18-1-4/cross.out.map")
    a.check()
