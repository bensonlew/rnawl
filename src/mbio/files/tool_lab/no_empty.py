# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class NoEmptyFile(File):
    '''
    非空文件
    '''
    def __init__(self):
        super(NoEmptyFile, self).__init__()

    def check(self):
        super(NoEmptyFile, self).check()
        if self.get_size() < 0:
            raise FileError('{} 文件为空'.format(self.path))
        return True
