# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.iofile import File
# import os
from biocluster.core.exceptions import FileError


class Level2File(File):
    """
    定义go_annotation比对输出go2level.xls格式
    """

    def __init__(self):
        super(Level2File, self).__init__()

    def get_info(self):
        super(Level2File, self).get_info()
        go2level_info = self.get_table_info()
        self.set_property('header', go2level_info[1])
        self.set_property('count', go2level_info[0])

    def get_table_info(self):
        with open(self.path) as f:
            go_info = f.read().split('\n')
            header = go_info[0]
            if len(header.split('\t')) != 6:
                raise FileError(
                    '文件表头存在错误，默认表头：term_type term GO number percent sequence', code = "43702801")
            count = 0
            for record in go_info[1:]:
                if record != '':
                    count += 1
                    record_info = record.split('\t')
                    if len(record_info) != 6:
                        raise FileError('文件中存在格式不整齐的行', code = "43702802")
                    if record_info[0] not in ['biological_process', 'cellular_component', 'molecular_function']:
                        raise FileError(
                            '文件中存在错误的一级分类：%s', variables = (record_info[0]), code = "43702803")
        return count, header

    def check(self):
        if super(Level2File, self).check():
            self.get_info()
            return True

    def get_gene(self):
        """获取gene列表"""
        with open(self.path) as f:
            f.readline()
            gene_list = []
            for line in f:
                line = line.strip('\n').split('\t')
                genes = line[-1].split(';')
                for i in genes:
                    gene_list.append(i.split('(')[0])
        return gene_list
