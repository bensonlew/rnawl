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
    核查vcf
        1.文本是否存在
        2.是否为空文本
        3.存在#注释内容，无变异位点数据
    self.path；self.prop["path"]
    """
    def __init__(self):
        super(VcfFile, self).__init__()

    def is_exists(self):
        if not os.path.isfile(self.path) or not os.path.exists(self.path):
            raise FileError("ERROR_:不存在%s文件！", variables=(self.prop["path"]), code="44800701")

    def get_info(self):
        with open(self.prop["path"], "r") as f:
            num = 0
            lines = f.readline().strip()
            while lines:
                if not re.match('#', lines):
                    num = 1     # 筛查只有#开头的内容,不然
                    # print("vcf核查通过{}".format(self.prop["path"]))
                    break       # 存在非#开头行的话，break
                lines = f.readline().strip()
            if num == 0:
                raise FileError("ERROR_:vcf%s文件无数据", variables=(self.prop["path"]), code="44800702")

    def check(self):
        super(VcfFile, self).check()
        self.is_exists()
        self.get_info()
        return True     # 可有可无


if __name__ == "__main__":
    a = VcfFile()
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/vcf_filecheck_data/null.vcf")
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/vcf_filecheck_data/sample2_0data.vcf")
    # a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/vcf_filecheck_data/right.vcf")
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/vcf_filecheck_data/0.vcf")
    a.check()
