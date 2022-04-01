# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class LgFile(File):
    """
    hongdong@20180627
    """

    def __init__(self):
        super(LgFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("原始文件中不存在%s文件！", variables=(self.path), code="44800201")

    def check(self):
        if super(LgFile, self).check():
            self.is_exists()
            return True
