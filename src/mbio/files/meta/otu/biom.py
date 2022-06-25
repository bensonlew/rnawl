# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""Biom格式文件类"""

from biocluster.iofile import File
import subprocess
import re
from biocluster.config import Config
import os
from biocluster.core.exceptions import FileError


class BiomFile(File):
    """
    Biom文件格式类, 需安装biom工具软件
    """
    def __init__(self):
        """
        """
        super(BiomFile, self).__init__()
        self.biom_path = os.path.join(Config().SOFTWARE_DIR, "miniconda2/bin/")

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(BiomFile, self).get_info()
        info = self.get_biominfo()
        self.set_property("form", self.biom_valid_check())
        self.set_property("otu_num", info[1])
        self.set_property("sample_num", info[1])
        self.set_property("metadata", info[2])

    def check(self):
        """
        检测文件是否满足要求
        :return:
        """
        if super(BiomFile, self).check():
            self.get_info()
            if self.prop['form']:
                pass
            else:
                raise FileError(u"文件格式错误", code="42700501")
        return True

    def biom_valid_check(self):
        """
        验证文件是否为biom格式
        """
        valid = False
        try:
            subpro = subprocess.check_output(self.biom_path+"biom  validate-table -i " + self.prop['path'], shell=True)
            if re.search("The input file is a valid BIOM-formatted file.", subpro):
                valid = True
        except subprocess.CalledProcessError:
            pass
        return valid

    def get_biominfo(self):
        """
        获取biom文件信息
        """
        try:
            subpro = subprocess.check_output("biom summarize-table -i " + self.prop['path'], shell=True)
            result = subpro.split('\n')
            for line in result:
                if line.startswith("Num samples"):
                    sample_num = re.split(r':\s+', line)[1]
                if line.startswith("Num observations"):
                    otu_num = re.split(r':\s+', line)[1]
                if line.startswith("Total count"):
                    seq_num = re.split(r':\s+', line)[1]
                if line.lstrip().startswith("Observation Metadata Categories"):
                    metadata = re.split(r':\s+', line)[1]
            return (sample_num, otu_num, seq_num, metadata)
        except subprocess.CalledProcessError:
            raise FileError("biom summarize-table 运行出错！", code="42700502")

    def convert_to_otu_table(self, otutable_filepath):
        """
        转换为OTU table格式
        """
        # biom convert -i otu_taxa_table.biom -o otu_taxa_table.txt --header-key taxonomy  --table-type \"OTU table\" --to-tsv
        # biom convert -i otu_table.biom -o otu_table.txt  --table-type \"otu table\"  --to-tsv
        cmd = self.biom_path+"biom convert -i " + self.prop['path'] + " -o " + otutable_filepath + ' --table-type \"OTU table\" --to-tsv'
        if self.prop['metadata'] == "taxonomy":
            cmd += " --header-key taxonomy "
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            raise FileError("biom convert 运行出错！", code="42700503")
        return True
