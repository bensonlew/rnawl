# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class SeqSampleFile(File):
    """
    定义 序列——样本 的文件格式
    这里的名称是序列名，当qc模块的输入是一个fastq文件的时候，用于拆分出各个样本
    """
    def __init__(self):
        super(SeqSampleFile, self).__init__()
        self.col = 0
        self.repeat_name = False

    def get_info(self):
        """
        获取文件属性
        """
        super(SeqSampleFile, self).get_info()

    def get_full_info(self):
        """
        因为这个文件较大，所以需要投递出去以后解析属性
        """
        (info, count) = self.get_file_info()
        self.set_property("sample_number", len(info))
        self.set_property("seq_number", count)
        self.set_property("sample_names", info.keys())

    def get_file_info(self):
        """
        获取seq_sample文件的信息
        """
        sample = dict()
        count = 0
        with open(self.prop['path'], 'r') as f:
            for line in f:
                count += 1
                line = line.rstrip('\n')
                line = re.split('\t', line)
                self.col = len(line)
                if line[1] not in sample.keys():
                    sample[line[1]] = 1
        return sample, count

    def full_check(self):
        if self.prop["sample_number"] == 0:
            raise FileError('应该至少包含一个样本')
        if self.col != 2:
            raise FileError('这个文件的列数为2')
        return True

    def check(self):
        if super(SeqSampleFile, self).check():
            self.get_info()
            return True
