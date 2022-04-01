# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class GoListFile(File):
    def __init__(self):
        super(GoListFile,self).__init__()

    def get_info(self):
        super(GoListFile,self).get_info()
        lines = open(self.prop['path']).readlines()
        for line in lines:
            items = line.strip().split('\t')
            if len(items) > 1:
                for i in items[1].split(';'):
                    if not i.startswith('GO'):
                        raise FileError('There is an incorrect comment -> {}'.format(i))

    def check(self):
        if super(GoListFile,self).check():
            self.get_info()
            return True