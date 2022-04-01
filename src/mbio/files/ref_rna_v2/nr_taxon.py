# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.iofile import File
# from biocluster.core.exceptions import FileError


class NrTaxonFile(File):
    """
    """
    def __init__(self):
        super(NrTaxonFile, self).__init__()

    def check(self):
        if super(NrTaxonFile, self).check():
            return True
