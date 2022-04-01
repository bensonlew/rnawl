# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180629
# src.mbio.files

from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os
import re


class VcfFile(File):
    """
    核查vcf文件是否存在
    """
    def __init__(self):
        super(VcfFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("ERROR_:不存在%s文件！", variables=(self.prop["path"]), code="44800701")

    def get_info(self):
        """
        不能为空文本
        存在#注释内容，无变异位点数据
        """
        with open(self.prop["path"], "r") as f:
            num = 0
            lines = f.readline().strip()
            while lines:
                if not re.match('#', lines):
                    num = 1     # 筛查只有#开头的内容,不然
                    break       # 存在非#开头行的话，break
                lines = f.readline().strip()
            if num == 0:
                raise FileError("ERROR_:vcf%s文件无数据", variables=(self.prop["path"]), code="44800702")

    def check(self):
        super(VcfFile, self).check()
        self.is_exists()


if __name__ == "__main__":
    a = VcfFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_evolution/file/pop.recode.vcf")
    a.check()
