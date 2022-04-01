# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class BedFile(File):
    """
    Binbinzhao@20180823
    """

    def __init__(self):
        super(BedFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在{}文件！".format(self.path))

    def check(self):
        if super(BedFile, self).check():
            self.is_exists()
            return True