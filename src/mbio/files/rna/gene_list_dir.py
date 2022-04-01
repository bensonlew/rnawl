# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError


class GeneListDirFile(Directory):
    """
    差异基因文件夹格式
    """
    def __init__(self):
        super(GeneListDirFile, self).__init__()

    def check(self):
        if super(GeneListDirFile, self).check():
            return True

    def get_info(self):
        super(GeneListDirFile, self).get_info()
        files = os.listdir(self.prop['path'])
        use_files = list()
        for f in files:
            if os.path.getsize(self.prop['path'] + '/' + f):
                use_files.append(self.prop['path'] + '/' + f)
        return use_files
