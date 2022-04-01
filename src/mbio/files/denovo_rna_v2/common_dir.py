# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
__author__ = "liubinxu"


class CommonDirFile(Directory):
    def __init__(self):
        super(CommonDirFile, self).__init__()

    def check(self):
        super(CommonDirFile, self).check()
