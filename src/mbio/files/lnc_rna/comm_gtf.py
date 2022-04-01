#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/14 15:46
@file    : comm_gtf.py
"""

from biocluster.iofile import File


class CommGtfFile(File):
    def __init__(self):
        super(CommGtfFile, self).__init__()
        self.__file_handler = None
        self.__file_type = None
        self.__line_parser = None

    def __init_handler(self):
        if self.__file_type is None:
            file_type = self.path.rsplit('.', 1)[-1].lower()
            self.set_type(file_type=file_type)
        self.__file_handler = self.get_reader()
        # self.__file_handler = open(r"C:\Users\zhipeng.zhao\Downloads\GCF_000001405.38_GRCh38.p12_genomic.gff")
        seek = 0
        while True:
            line = self.__file_handler.readline()
            if not line.startswith('#') and line != '\n':
                self.__file_handler.seek(seek)
                break
            seek = self.__file_handler.tell()

    def set_type(self, file_type='gtf'):
        """set the type of file to assign line parser to self.__line_parser

        :param file_type: gtf or gff
        :return:
        """
        if file_type == 'gtf':
            self.__line_parser = self.__get_parser(items_sep=' ')
        elif file_type == 'gff' or file_type == 'gff3':
            self.__line_parser = self.__get_parser(items_sep='=')
        else:
            raise Exception('file type set error, it should be gtf or gff')

    def __get_parser(self, items_sep=' '):
        def _parser(line):
            line_splits = line.strip('; \r\n').split('\t')
            line_splits[3] = int(line_splits[3])
            line_splits[4] = int(line_splits[4])
            if len(line_splits[8].split('; ')) >= 2:
                line_splits[8] = {
                    k.strip(): v.strip(' "')
                    for k, v in (
                        item.strip().split(items_sep, 1) for item in line_splits[8].split('; ')
                    )
                }
            else:
                line_splits[8] = {
                    k.strip(): v.strip(' "')
                    for k, v in (
                        item.strip().split(items_sep, 1) for item in line_splits[8].split(';')
                    )
                }
            return line_splits
        return _parser

    def set_path(self, path):
        super(CommGtfFile, self).set_path(path)
        return self

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
        if line.startswith('#') or line == '\n':
            while True:
                line = self.__file_handler.readline()
                if line.startswith('#') or line == '\n':
                    continue
                else:
                    break
        if not line:
            raise StopIteration

        line_splits = self.__line_parser(line)
        return line_splits, line

    def next(self):
        """
        python2 iterator iter method
        :return:
        """
        return self.__next__()


if __name__ == '__main__':
    with CommGtfFile() as lnc_file:
        lnc_file.set_type('gff')
        for line in lnc_file:
            print(line[0])
            print(line[1])
