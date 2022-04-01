# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class ExpressMatrixFile(File):
    """
    """

    def __init__(self):
        super(ExpressMatrixFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(ExpressMatrixFile, self).get_info()
        samples, genes = self.get_matrix_info()
        self.set_property('sample', samples)
        self.set_property('gene', genes)

    def get_matrix_info(self):
        """
        获取并返回表达量矩阵的信息
        :return:
        """
        with open(self.prop['path']) as tempfile:
            lines = tempfile.readlines()
            genes = []
            if not lines:
                raise FileError('表达量矩阵为空')
            for line in lines[1:]:
                genes.append(line.strip('\n').split('\t')[0])
            samples = lines[0].strip('\n').split('\t')[1:]
        return samples, genes

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        if super(ExpressMatrixFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            self.get_info()
            if self.prop['gene'] == 0:
                raise FileError('表达量矩阵至少有一个基因')

    def get_list(self, output):
        with open(self.prop['path'], 'rb') as r, open(output, 'wb') as w:
            i = 0
            for f in r:
                i += 1
                if i == 1:
                    pass
                else:
                    w.write(f.split('\t')[0] + '\n')
        return output
