# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import os
import re
from biocluster.core.exceptions import FileError
from biocluster.iofile import Directory
from mbio.files.meta.otu.otu_table import OtuTableFile
from mbio.files.meta.otu.biom import BiomFile


class TaxSummaryDirFile(Directory):
    """
    定义tax_summary_dir文件夹格式
    """
    def __init__(self):
        super(TaxSummaryDirFile, self).__init__()
        self.biom = 0
        self.otu_table = 0

    def get_info(self):
        """
        获取文件夹属性
        """
        if 'path' in self.prop.keys() and os.path.isdir(self.prop['path']):
            self.get_file_number()
            self.set_property('biom_number', self.biom)
            self.set_property('otu_number', self.otu_table)

    def get_file_number(self):
        """
        获取文件夹下的biom文件和otu_table文件的数目
        :return:文件数目
        """
        filelist = os.listdir(self.prop['path'])
        for file_ in filelist:
            file_ = os.path.join(self.prop['path'], file_)
            otu = OtuTableFile()
            try:
                otu.set_path(file_)
                self.otu_table += 1
            except FileError:
                pass
            biom = BiomFile()
            try:
                biom.set_path(file_)
                self.biom += 1
            except FileError:
                pass
        return (self.biom, self.otu_table)

    def get_table(self, level, full_path=False):
        """
        获取OTU表
        :param level: 输入的等级
        :param full_path: 是要full.xls表还是stat.xls表
        """
        list_ = os.listdir(self.prop['path'])
        for file_ in list_:
            if full_path:
                pattern = r"otu_taxon_" + level + r"\.full\.xls"
                if re.search(pattern, file_, re.IGNORECASE):
                    return os.path.join(self.prop['path'], file_)
            else:
                pattern = r"otu_taxon_" + level + r"\.xls"
                if re.search(pattern, file_, re.IGNORECASE):
                    return os.path.join(self.prop['path'], file_)
        raise ValueError("未找到文件，输入的level为:{}，正确的level应该为['Domain', 'Kingdom', 'Phylum', 'Order', 'Family', 'Genus','Class' , 'Species', 'otu']中的一个".format(level))

    def get_biom(self, level):
        """
        获取biom表
        :param level: 输入的等级
        """
        list_ = os.listdir(self.prop['path'])
        for file_ in list_:
            pattern = r"otu_taxon_" + level + r"\.biom"
            if re.search(pattern, file_, re.IGNORECASE):
                return os.path.join(self.prop['path'], file_)
        raise ValueError("未找到文件，检查输入的level是否正确")

    def check(self):
        """
        检测文件夹是否满足要求，不满足时触发FileError异常
        """
        # self.get_info()
        if super(TaxSummaryDirFile, self).check():
            self.get_info()
            if not os.path.isdir(self.prop['path']):
                raise FileError("不是一个文件夹", code="42701001")
                return False
            return True
            # if self.biom % 8 != 0:
            #    raise FileError("文件格式不正确")

if __name__ == "__main__":
    a = TaxSummaryDirFile()
    a.set_path("tax_dir")
    a.check()
    print a.prop
