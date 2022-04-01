# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""OTU seqids格式文件类"""

from biocluster.iofile import File
import subprocess


class OtuSeqidsFile(File):
    """
    OtuSeqids
    """
    def __init__(self):
        """
        """
        super(OtuSeqidsFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        :return:
        """
        super(OtuSeqidsFile, self).get_info()
        self.set_property("otu_num", self.get_otunum())

    def get_otunum(self):
        """
        获取otu数目
        """
        try:
            subpro = subprocess.check_output("wc -l " + self.prop['path'], shell=True)
            return int(subpro.split(' ')[0])
        except subprocess.CalledProcessError:
            raise Exception("wc -l 运行出错！")
