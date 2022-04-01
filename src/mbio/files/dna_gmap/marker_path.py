# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.07.02

import os
import re
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class MarkerPathFile(Directory):
    """
    文件夹中必须要有map和marker文件。
    """
    def __init__(self):
        super(MarkerPathFile, self).__init__()

    def dir_check(self):
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            n = 0
            for m in os.listdir(self.prop['path']):
                if re.match(r".*\.pri\.map$", m):
                    n += 1
                if re.match(r".*\.pri\.marker$", m):
                    n += 1
            if n == 0:
                raise FileError("文件夹路径中没有map和marker文件!", code="44800901")

    def check(self):
        if super(MarkerPathFile, self).check():
            self.dir_check()
            return True

if __name__ == "__main__":
    a = MarkerPathFile()
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/gmap/05.map-cycle1-18-1-4")
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/gmap")
    a.check()
