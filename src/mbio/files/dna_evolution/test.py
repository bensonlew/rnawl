# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.03

import os
from biocluster.iofile import File
from biocluster.core.exceptions import FileError


class TestFile(File):
    """
    测试
    """
    def __init__(self):
        super(TestFile, self).__init__()

    def check(self):
        if super(TestFile, self).check():
            return True


if __name__ == "__main__":
    a = TestFile()
    a.set_path("/mnt/ilustre/users/sanger-dev/workspace/20180920/Single_sweep_region_single/SweepRegionSingle/LdRegion/output/1-2.2.pi_tajimaD_fst.select.region")
    a.check()
