#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/1 16:49
@file    : lnc_gtf.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
from collections import OrderedDict

from biocluster.iofile import File


class LncGtfFile(File):
    def __init__(self, file_path=None):
        super(LncGtfFile, self).__init__()
        if file_path is not None:
            self.set_path(file_path)
        self.__file_handler = None

    def __init_handler(self):
        # self.__file_handler = open(r"C:\Users\zhipeng.zhao\Desktop\old_transcripts.gtf")
        self.__file_handler = self.get_reader()
        seek = 0
        while True:
            line = self.__file_handler.readline()
            if not line.startswith('#') and line != '\n':
                self.__file_handler.seek(seek)
                break
            seek = self.__file_handler.tell()

    def __enter__(self):
        if self.__file_handler is None:
            self.__init_handler()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()

    def __iter__(self):
        if self.__file_handler is None:
            self.__init_handler()
        return self

    def __next__(self):
        """
        python3 iterator iter method
        :return:
        """
        if self.__file_handler is None:
            self.__init_handler()
        line = self.__file_handler.readline()
        if not line:
            raise StopIteration
        line_splits = line.strip().split('\t')
        line_splits[3] = int(line_splits[3])
        line_splits[4] = int(line_splits[4])
        if len(line_splits[8].split('; ')) >= 2:
            line_splits[8] = OrderedDict(
                (k.strip(), v.strip(' "'))
                for k, v in (
                item.strip().split(' ', 1) for item in line_splits[8].strip('; ').split('; ')))
        else:
            line_splits[8] = OrderedDict(
                (k.strip(), v.strip(' "'))
                for k, v in (
                item.strip().split(' ', 1) for item in line_splits[8].strip('; ').split(';')))

        return line_splits, line

    def next(self):
        """
        python2 iterator iter method
        :return:
        """
        return self.__next__()


if __name__ == '__main__':
    for line in LncGtfFile():
        print(line[0])
        print(line[1])
