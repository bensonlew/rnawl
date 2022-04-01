# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from mbio.files.small_rna.common import CommonFile
import ConfigParser

__author__ = "liubinxu"

class IniFile(CommonFile):
    def __init__(self):
        super(IniFile, self).__init__()

    def check(self):
        super(IniFile, self).check()
