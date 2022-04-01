# -*- coding: utf-8 -*-
# __author__ = gaohao
# last_modify：2019.01.15


from biocluster.iofile import Directory
import os,re
from biocluster.core.exceptions import FileError

class InputDirFile(Directory):
    """
    输入文件转化成PAGP的输入文件
    """
    def __init__(self):
        super(InputDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(InputDirFile, self).get_info()
        file_list = os.listdir(self.prop['path'])
        for file in  file_list:
            pass
            # if not file.endswith((".ptt", ".ffn", ".faa", ".fna",".fasta",'.txt')):
            #     raise FileError('%s文件后缀格式不正确，请检查！', variables=(file))


    def check(self):
        if super(InputDirFile, self).check():
            self.get_info()
            return True

if __name__ == "__main__":
    dir = InputDirFile()
    dir.check()