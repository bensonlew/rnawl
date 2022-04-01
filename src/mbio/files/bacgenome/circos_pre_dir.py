# -*- coding: utf-8 -*-
# __author__ = zouxuan
# last_modify：2018.04.09


from biocluster.iofile import Directory
import os
from biocluster.core.exceptions import FileError

'''
circos准备文件夹的检查
'''


class CircosPreDirFile(Directory):
    """
    """

    def __init__(self):
        super(CircosPreDirFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(CircosPreDirFile, self).get_info()
        file_list = os.listdir(self.prop['path'])
        for file in ["karyotype.txt", "sense_strand_cog.txt", "temp.txt", "antisense_strand_cog.txt",
                     "positive_gc_count.txt", "negative_gc_count.txt", "positive_gc_skew.txt", "negative_gc_skew.txt"]:
            if file not in file_list:
                raise FileError('%s文件不存在，请检查', variables=(file), code="41400301")
            self.set_property("file_list", file_list)

    def check(self):
        if super(CircosPreDirFile, self).check():
            self.get_info()
            return True


if __name__ == "__main__":
    dir = CircosPreDirFile()
