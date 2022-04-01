# -*- coding: utf-8 -*-
"""
@time    : 2018/10/25 11:04
@file    : diff_gene_list.py
@author  : zhipeng.zhao
@contact: 757049042@qq.com
"""
from biocluster.iofile import File

class DiffGeneListFile(File):
    """
    差异基因列表文件
    """
    def __init__(self):
        super(DiffGeneListFile, self).__init__()

    def get_info(self):
        """
        文件属性获取
        """
        super(DiffGeneListFile, self).get_info()
        gene_list, gene_num = self.gene_list
        self.set_property('gene_list', gene_list)
        self.set_property('gene_num', gene_num)

    def check(self):
        if super(DiffGeneListFile, self).check():
            self.get_info()
            return True

    @property
    def gene_list(self):
        gene_list = []
        with self.get_reader() as infile:
            gene_list = [
                i.strip() for i in infile.read().strip('\n').split("\n")
            ]
        length = len(gene_list)
        return gene_list, length
