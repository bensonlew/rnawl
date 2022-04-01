# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.iofile import File
from biocluster.core.exceptions import FileError
import os


class QzvFile(File):
    """
    定义metaasv的qzv格式的文件格式
    """
    def __init__(self):
        super(QzvFile, self).__init__()

    def check_info(self):
        """
        获取判断文件名称
        """
        file_name = os.path.basename(self.prop['path'])
        if not file_name.endswith(".qzv"):
            raise FileError('文件不是以qzv结尾，请检查文件类型！')
        return True

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        if super(QzvFile, self).check():
            self.check_info()
            return True
