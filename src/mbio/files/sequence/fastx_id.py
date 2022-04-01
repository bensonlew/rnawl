# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.iofile import File
import re
from biocluster.core.exceptions import FileError


class FastxIdFile(File):
    """
    fasta和fastq的id文件
    """

    def __init__(self):
        super(FastxIdFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(FastxIdFile, self).get_info()
        fastx_id_info = self.get_fastx_id_info()
        self.set_property('fastx_id_list', fastx_id_info)

    def get_fastx_id_info(self):
        """
        获取并返回id列表
        :return:
        """
        fastx_id_list = []
        try:
            with open(self.prop['path']) as idfile:
                for fastx_id in idfile.readlines():
                    if re.match('[\t\n]+', fastx_id):
                        pass
                    else:
                        fastx_id_list.append(fastx_id)
        except IOError:
            raise FileError('无法打开文件')
        return fastx_id_list

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        :return:
        """
        self.get_info()
        if super(FastxIdFile, self).check():
            raise FileError('文件格式不正确')
