# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class SgrnaRegionFile(File):
    """
    zhaobinbin@20190320
    """

    def __init__(self):
        super(SgrnaRegionFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在{}文件！".format(self.path))

    def check(self):
        if super(SgrnaRegionFile, self).check():
            self.is_exists()
            return True