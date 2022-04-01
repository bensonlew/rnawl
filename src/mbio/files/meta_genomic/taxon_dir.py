# -*- coding: utf-8 -*-
import re
import os
import glob
from collections import defaultdict
from biocluster.config import Config
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile


levels = ['D', 'K', 'P', 'C', 'O', 'F', 'G', 'S']


class TaxonDirFile(Directory):
    """
    定义宏基因sec注释结果文件夹的格式
    """
    def __init__(self):
        super(TaxonDirFile, self).__init__()
        self.taxon_list = []
        self.level_list = defaultdict(list)

    def get_taxon_result(self):
        self.taxon_list = glob.glob(os.path.join(self.path, "*.taxon.xls"))
        self.set_property("taxon_list", self.taxon_list)
        print(self.taxon_list)

    def get_level_list(self):
        for index, level in enumerate(levels):
            post_fix = "*_{}.xls".format(level)
            self.level_list[str(index + 1)] = glob.glob(os.path.join(self.path, post_fix))
        self.set_property("level_list", self.level_list)
        print(self.level_list)

    def get_level(self, level_id):
        return self.level_list[str(level_id)][0]

    def check_info(self):
        """
        获取判断文件夹内容
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            file_name = os.listdir(self.prop["path"])
            if len(file_name) == 0:
                raise FileError("文件夹为空，请设置正确的文件夹路径！", code="45301001")
            self.get_taxon_result()
            self.get_level_list()
        else:
            raise FileError("文件夹路径不正确，请设置正确的文件夹路径!", code="45301002")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        :return:
        """
        super(TaxonDirFile, self).check()
        self.check_info()
        return True
