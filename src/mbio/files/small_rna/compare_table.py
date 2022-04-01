# -*- coding: utf-8 -*-
# __author__ = 'gdq'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError

class CompareTableFile(File):
    """
    cmp_info: path of cmp info, file with only two columns.
              Header line starts with '#', or Just give no header line.
              lines separated with tab/blank.
              group/sample name should be a word.
               -----------------
               #ctrl    test
               group1   group2
               group2   group3
               s5       s6
               -----------------
    """
    def __init__(self):
        super(CompareTableFile, self).__init__()

    def get_info(self):
        super(CompareTableFile, self).get_info()
        cmp_list = self.parse_file()
        self.set_property('cmp_list', cmp_list)

    def parse_file(self):
        # comparison info -> list. [(ctrl, test), ...]
        with open(self.prop['path']) as f:
            cmp_list = list()
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                tmp_ctrl, tmp_test = line.strip().split()
                cmp_list.append((tmp_ctrl, tmp_test))
        cmp_list = sorted(list(set(cmp_list)))
        if not cmp_list:
            raise FileError('比较信息的内容为空', code = "42000601")
        return cmp_list

    def check(self):
        super(CompareTableFile, self).check()
        self.get_info()
