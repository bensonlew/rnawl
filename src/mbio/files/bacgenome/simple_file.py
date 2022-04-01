# -*- coding: utf-8 -*-
import os
import pandas as pd
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
__author__ = "guhaidong"


class SimpleFileFile(File):
    def __init__(self):
        super(SimpleFileFile, self).__init__()

    def check(self):
        super(SimpleFileFile, self).check()
