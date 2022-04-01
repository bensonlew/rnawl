# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
__author__ = "liubinxu"


class CommonFile(File):
    def __init__(self):
        super(CommonFile, self).__init__()

    def check(self):
        super(CommonFile, self).check()
