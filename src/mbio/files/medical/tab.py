# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# creat at 20171111

from biocluster.iofile import File
import re
from biocluster.core.exceptions import FileError


class TabFile(File):
    '''
    tab或逗号分割的表格文件, 有表头
    '''
    def __init__(self):
        super(TabFile, self).__init__()
        self.rows = list()
        self.colums = list()

    def get_info(self):
        '''
        获取文件属性
        '''
        super(TabFile, self).get_info()
        with open(self.prop['path'], "r") as f:
            for line in f.readlines():
                cols = re.split(r'[,\t]', line.strip())
                self.rows.append(cols)
        self.colums = zip(*self.rows)

    def check(self):
        '''
        检查文件格式
        '''
        if super(TabFile, self).check():
            self.get_info()
        else:
            raise FileError("文件格式错误")
        return True

    def get_colums(self, ids):
        '''
        按列的index获取多列的list
        '''
        return [self.colums[i] for i in ids]

    def get_rows(self, ids):
        '''
        按列的index获取多列的list
        '''
        return [self.rows[i] for i in ids]
