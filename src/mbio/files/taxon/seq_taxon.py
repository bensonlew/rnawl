# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""Taxon格式文件类"""

from biocluster.iofile import File
import re
from biocluster.core.exceptions import FileError
from collections import defaultdict


class SeqTaxonFile(File):
    """
    Taxon文件格式类
    """
    def __init__(self):
        """
        """
        super(SeqTaxonFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(SeqTaxonFile, self).get_info()
        info = self.get_fileinfo()
        self.set_property("form", info[0])
        self.set_property("seq_num", info[1])

    def check(self):
        """
        检测文件是否满足要求
        :return:
        """
        if super(SeqTaxonFile, self).check():
            self.get_info()
            if self.prop['form']:
                pass
            else:
                raise FileError(u"文件格式错误", code="44200101")
        return True

    def get_fileinfo(self):
        """
        获取物种分类文件信息
        """
        form, seq_num = True, 0
        with open(self.prop['path'], 'r') as f:
            while 1:
                line = f.readline().rstrip()
                if not line:
                    break
                taxline = re.split(r'\t', line)[1]
                taxs = re.split(r';\s*', taxline)
                # if len(taxs) != 8:
                #     raise FileError("%s序列注释信息不为8层[dkpcofgs]，请检查输入的taxon文件！", variables=(line[0]), code="44200103")
                for t in taxs:
                    if re.match(r'[dkpcofgs]\_\_\S+', t):
                        pass
                    else:
                        form = False
                        break
                seq_num += 1
                if not line:
                    break
        return (form, seq_num)

    def get_all_name(self):
        my_name = defaultdict(int)
        with open(self.prop["path"], 'rb') as r:
            for line in r:
                line = re.split('\t', line)
                my_name[line[0]] += 1
        dup_list = list()
        for k in my_name.iterkeys():
            if my_name[k] > 1:
                dup_list.append(k)
        if len(dup_list) > 0:
            str_ = "; ".join(dup_list)
            raise FileError("序列名:%s在输入的tax文件里面重复", variables=(str_), code="44200102")
        return my_name
