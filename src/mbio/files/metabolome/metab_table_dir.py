# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re,os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
from collections import defaultdict
from metab_desc import MetabDescFile
from metab_abun import MetabAbunFile


class MetabTableDirFile(Directory):
    """
    定义metab_table格式文件
    """
    def __init__(self):
        super(MetabTableDirFile, self).__init__()
        self.file_list = ['metab_abund.txt', 'metab_desc.txt']

    def get_info(self):
        """
        获取文件属性
        """
        super(MetabTableDirFile, self).get_info()
        if not os.path.isdir(self.prop['path']):
            raise FileError("文件夹路径不正确，%s" , variables=(self.prop['path']), code="44700901")
        for file in self.file_list:
            file_rout = os.path.join(self.prop['path'], file)
            if not os.path.isfile(file_rout):
                raise FileError("文件夹中不存在文件：%s,请检查" , variables=(file_rout), code="44700902")
        info = self.get_file_info()
        self.set_property("metab_abun", info[0])
        self.set_property("metab_desc", info[1])
        # self.set_property("sample_number", len(info[0]))
        # self.set_property("sample", info[0])
        # self.set_property("group_scheme", info[1])
        # self.set_property("is_empty", info[2])

    def get_file_info(self):
        """
        获取metab_table文件的信息
        """
        abun_obj = MetabAbunFile()
        abun_obj.set_path(os.path.join(self.prop['path'], "metab_abund.txt"))
        desc_obj = MetabDescFile()
        desc_obj.set_path(os.path.join(self.prop['path'], "metab_desc.txt"))
        return abun_obj,desc_obj

    def format_check(self):
        self.check_abun()
        self.check_desc()

    def check_abun(self):
        if self.prop.has_key('metab_abun'):
            self.prop['metab_abun'].check()
        else:
            raise FileError("没有metab_abun属性，先运行get_info方法", code="44700903")

    def check_desc(self):
        if self.prop.has_key('metab_desc'):
            self.prop['metab_desc'].check()
        else:
            raise FileError("没有metab_desc属性，先运行get_info方法", code="44700904")

    def check(self):
        if super(MetabTableDirFile, self).check():
            self.get_info()
            self.format_check()

if __name__ == "__main__":
    g = MetabTableDirFile()
    g.set_path("example.group")
    g.get_info()
    g.sub_group("example.group.sub", ["g1", "g3", "g4"])
