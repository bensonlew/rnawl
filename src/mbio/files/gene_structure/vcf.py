# -*- coding: utf-8 -*-
# __author__ = 'chenyy'
import os
import re
import gzip
import subprocess
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.config import Config

class VcfFile(File):
    """
    定义vcf文件
    """
    def __init__(self):
        super(VcfFile, self).__init__()
          
    def check(self):
        """
        检查文件是否满足要求，不满足触发FileError异常
        """
        if super(VcfFile, self).check():
            return True
        else:
            raise FileError("文件格式错误！")
    """
        
     def get_info(self):
       
        super(VcfFile, self).get_info()
        self.get_vcf_info()


        
    def get_vcf_info(self):
        with open(self.prop["path"],"r") as v:
            lines = v.readlines()
            
    """