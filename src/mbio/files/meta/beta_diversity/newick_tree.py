# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.iofile import File
import re
# from Bio import Phylo
# import subprocess
# from biocluster.config import Config
# import os
from biocluster.core.exceptions import FileError


class NewickTreeFile(File):
    """
    """

    def __init__(self):
        super(NewickTreeFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(NewickTreeFile, self).get_info()
        newickinfo = self.get_newick_info()
        self.set_property('sample', newickinfo[:: -1])
        self.set_property('count', len(newickinfo))

    def get_newick_info(self):
        """
        获取并返回树文件信息
        :return:
        """
        tempfile = open(self.prop['path'])
        lines = tempfile.readlines()
        if not lines:
            raise FileError('树文件为空', code="42700301")
        tree = lines[0].rstrip()
        raw_samp = re.findall(r'([(,]([\[\]\.\;\'\"\ \[\]0-9a-zA-Z_-]+?):[0-9])', tree)
        samp = [i[1] for i in raw_samp]
        return samp

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(NewickTreeFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            tempfile = open(self.prop['path'])
            tree = tempfile.readlines()
            if len(tree) == 1:
                pass
            else:
                raise FileError('文件中存在多个newick树', code="42700302")
                # 可以保留第一个数，删除其他树，继续操作
            tree = tree[0].rstrip()
            if tree.count('(') == tree.count(')'):
                pass
            else:
                raise FileError('树文件格式错误', code="42700303")
            if len(re.findall(r':[\.0-9]+', tree)) == tree.count(':'):
                pass
            else:
                raise FileError('程序只接受带有分支距离的树文件，或者文件中距离表示错误', code="42700304")
            if tree[-1] == ';':
                pass
            else:
                raise FileError('文件结尾不是分号‘;’', code="42700305")
        return True

    # def terminals_rank(self):
    #     """
    #     获取树的枝叶排序列表
    #     """
    #     self.check()
    #     tree = Phylo.read(self.path, 'newick')
    #     terminals = [i.name for i in tree.get_terminals()]
    #     return terminals[:: -1]
