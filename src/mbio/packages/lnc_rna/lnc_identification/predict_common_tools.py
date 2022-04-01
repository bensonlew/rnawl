#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 9:55
@file    : predict_common_tools.py
"""


class FastaParser(object):
    def __init__(self, file, record_head_parser=None):
        self.__file = file
        self.__file_handler = None
        self.__head_parser = record_head_parser if callable(record_head_parser) else lambda x: x.strip()

    def __iter__(self):
        if self.__file_handler is None:
            self.__file_handler = open(self.__file)
        flag = 0
        record_dic = {}
        fa_seq_list = []
        for line in self.__file_handler:
            if line.startswith('>'):
                if flag == 0:
                    record_dic['head'] = self.__head_parser(line)
                    flag += 1
                    continue
                record_dic['fa_seq'] = ''.join(fa_seq_list)
                yield record_dic
                fa_seq_list = []
                record_dic = {'head': self.__head_parser(line)}
            else:
                fa_seq_list.append(line)
        record_dic['fa_seq'] = ''.join(fa_seq_list)
        yield record_dic

    def __enter__(self):
        self.__file_handler = open(self.__file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()


class GFFParser(object):
    def __init__(self, file):
        self.__file = file
        self.__file_handler = None

    def __iter__(self):
        if self.__file_handler is None:
            self.__file_handler = open(self.__file)
        for line in self.__file_handler:
            if line.startswith('#'):
                continue
            line_splits  = line.strip().split('\t')
            line_splits[3] = int(line_splits[3])
            line_splits[4] = int(line_splits[4])
            line_splits[8] = {
                k.strip(): v.strip(' "')
                for k, v in (
                    item.strip().split('=', 1) for item in line_splits[8].strip('; ').split(';')
                )
            }
            yield line_splits, line

    def __enter__(self):
        self.__file_handler = open(self.__file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()


class GTFParser(object):
    def __init__(self, file):
        """

        :param file:
        :param return_record: True --> return transcript record, False --> return line info
        """
        self.__file = file
        self.__file_handler = None

    def __iter__(self):
        if self.__file_handler is None:
            self.__file_handler = open(self.__file)
        for line in self.__file_handler:
            if line.startswith('#'):
                continue
            line_splits  = line.strip().split('\t')
            line_splits[3] = int(line_splits[3])
            line_splits[4] = int(line_splits[4])
            line_splits[8] = {
                k.strip(): v.strip(' "')
                for k, v in (
                    item.strip().split(' ', 1) for item in line_splits[8].strip('; ').split(';')
                )
            }
            yield line_splits, line

    def __enter__(self):
        self.__file_handler = open(self.__file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()


class CommandFactory(object):
    def __init__(self, exe_path):
        """
        :param exe_path: python cpat.py
        """
        self.__params_list = [exe_path]
        self.__cmd = None

    def add_param(self, param, param_value, param_desc=None):
        self.__params_list.extend((param, param_value))

    def to_string(self):
        if self.__cmd is None:
            self.__cmd = ' '.join(str(i) for i in self.__params_list)
        return self.__cmd

    def __str__(self):
        return self.to_string()


if __name__ == '__main__':
    with FastaParser(r"C:\Users\zhipeng.zhao\Desktop\test_lnc_predict.fa") as fa_handler:
        for dic in fa_handler:
            print(dic)
