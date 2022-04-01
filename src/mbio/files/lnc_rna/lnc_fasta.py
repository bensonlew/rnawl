#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/1 15:04
@file    : lnc_fasta.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""

from biocluster.iofile import File
from Bio import SeqIO


class LncFastaFile(File):
    def __init__(self, file_path_=None):
        super(LncFastaFile, self).__init__()
        self.__file_handler = None
        self.seq_list = []
        self.seq_head = None
        if file_path_ is not None:
            self.set_path(file_path_)

    def __init_handler(self):
        # self.__file_handler = open(r'C:/Users/zhipeng.zhao/Desktop/test_lnc_predict.fa')
        self.__file_handler = self.get_reader()
        self.seqio = SeqIO.parse(self.path, "fasta")

        '''
        for line in self.__file_handler:
            if line == '\n' or line.startswith('#'):
                continue
            elif line.startswith('>'):
                self.seq_head = line
                break
        '''

    def __enter__(self):
        if self.__file_handler is None:
            self.__init_handler()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()

    def __iter__(self):
        if self.__file_handler is None:
            self.__init_handler()
        for seq in self.seqio:
            yield(seq.id, seq.seq)

    '''
    def __next__(self):
        """
        python3 iterator iter method
        :return:
        """
        # for seq in self.seqio:
           #  yield (seq.id, seq.seq)
        if self.seq_head is None:
            raise StopIteration
        if self.__file_handler is None:
            # print self.seq_head
            self.__init_handler()

        while True:
            # print line
            line = self.__file_handler.readline()
            if line.startswith('>') or (not line and len(self.seq_list) > 0):
                temp_list = self.seq_list
                temp_head = self.seq_head
                self.seq_list = []
                self.seq_head = line
                return temp_head, ''.join(temp_list)
            elif not line:
                return
                # raise StopIteration
            self.seq_list.append(line)
        '''
    '''
    def next(self):
        """
        python2 iterator iter method
        :return:
        """
        return self.__next__()
    '''


if __name__ == '__main__':
    # handler = LncFastaFile()
    # handler.set_path(r'C:/Users/zhipeng.zhao/Desktop/test_lnc_predict.fa')
    with LncFastaFile() as f:
        for line in f:
            print(line[0])
            print(line[1])
