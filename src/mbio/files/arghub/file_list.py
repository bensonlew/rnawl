# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.iofile import File
import os
from biocluster.core.exceptions import FileError


class FileListFile(File):
    """
    """
    def __init__(self):
        super(FileListFile, self).__init__()
        self.file_list = []

    def check(self):
        if super(FileListFile, self).check():
            err_path = []
            with open(self.path, 'r') as r:
                for l in r:
                    line = l.strip().split()
                    if os.path.isfile(line[-1]):
                        self.file_list.append(line[-1])
                    else:
                        err_path.append(line[-1])
            if err_path:
                raise FileError('文件list {}中，一下文件不存在: {}'.format(self.path, err_path))
        return True

