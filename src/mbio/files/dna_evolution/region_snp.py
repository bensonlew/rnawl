# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.03

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class RegionSnpFile(File):
    """
    测试
    """
    def __init__(self):
        super(RegionSnpFile, self).__init__()

    def check(self):
        if super(RegionSnpFile, self).check():
            is_null = self.get_info()
            self.set_property("is_null", is_null)

    def get_info(self):
        """
        检查文件是否是空的
        """
        with open(self.prop["path"], "r") as f:
            lines = f.readlines()
            if len(lines) == 1:
                return True
        return False


if __name__ == "__main__":
    a = RegionSnpFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_tsg_33071_SweepRegion_1220094354535412/SweepRegion/SweepRegionSingle/FstManhattan/output/Q1-Q3.fst_pi.detail.select")
    a.check()
