# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class QualFile(File):
    """
    定义Fasta文件的质量文件类型
    """
    def __init__(self):
        super(QualFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(QualFile, self).get_info()
        seqinfo = self.get_seq_info()
        self.set_property("seq_number", seqinfo)

    def check(self):
        """
        检测文件是否满足要求,发生错误时应该触发FileError异常
        """
        self.get_info()
        if super(QualFile, self).check():
            try:
                qualfile = open(self.prop['path'])
            except IOError:
                raise FileError('无法打开文件')
            line1 = qualfile.readline()
            if line1[0] != '>':
                raise FileError('文件格式错误')
            line2 = qualfile.readline().rstrip().split()
            for base in line2:
                try:
                    int(base)
                except ValueError:
                    raise FileError('错误的文件格式')
            qualfile.close()
        return True

    def get_seq_info(self):
        """
        获取qual信息
        """
        qualfile = open(self.prop['path'])
        num = 0
        for line in qualfile:
            if line[0] == '>':
                num += 1
        return num
