# -*- coding: utf-8 -*-
# __author__ = 'chenhongyu'

'''PT_NIPT 实验批次表、客户信息表、samples表'''

from biocluster.iofile import File
import subprocess
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import FileError


class SamplesRefFile(File):
    def __init__(self):
        super(SamplesRefFile, self).__init__()

    def check(self):
        if super(SamplesRefFile, self).check():
            return True
        else:
            raise FileError("文件未设置")











