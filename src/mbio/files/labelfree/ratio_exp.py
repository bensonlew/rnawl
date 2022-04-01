# -*- coding: utf-8 -*-
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from collections import OrderedDict
__author__ = 'shicaiping'


class RatioExpFile(File):

    def __init__(self):
        super(RatioExpFile, self).__init__()

    def get_info(self):
        super(RatioExpFile, self).get_info()
        samples = self.parse_file()
        self.set_property("sample_number", len(samples))
        self.set_property("samples", samples)

    def parse_file(self):
        sample_list = set()
        with open(self.prop['path']) as f:
            header_line = f.readline()
        tmp_list = header_line.strip().split()
        for s in tmp_list[1:]:
            sample_list.add(s)
        return sorted(sample_list)

    def check(self):
        super(RatioExpFile, self).check()
        self.get_info()
