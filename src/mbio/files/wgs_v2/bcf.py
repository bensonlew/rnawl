# -*- coding: utf-8 -*-

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class BcfFile(File):
    """
    wentian@20190220
    """

    def __init__(self):
        super(BcfFile, self).__init__()

    def check(self):
        if super(BcfFile, self).check():
            if not os.path.exists(self.prop["path"]):
                raise FileError("%s文件不存在" % self.prop["path"])
            return True
