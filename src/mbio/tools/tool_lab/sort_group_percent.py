# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import os
import re
import xlrd
import shutil

class SortGroupPercentAgent(Agent):
    """
    对应关系circos图用于生成组内占百分比表
    """
    def __init__(self, parent):
        super(SortGroupPercentAgent, self).__init__(parent)
        options = [
            {"name":"table_file","type":"infile","format":"small_rna.common"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("table_file"):
            raise OptionError("必须提供表格文件")
    
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SortGroupPercentAgent, self).end()

class SortGroupPercentTool(Tool):
    def __init__(self,config):
        super(SortGroupPercentTool, self).__init__(config)
        self._version = 'v.10'

    def run(self):
        """
        运行小工具
        """
        super(SortGroupPercentTool, self).run()
        group_table = self.option("table_file").prop["path"]
        out_table = os.path.join(self.output_dir,"group_percents_table.txt")
        with open(group_table,"r") as gt, open(out_table,"w") as ot:
            first_title = gt.readline()
            ot.write(first_title)
            while 1:
                line = gt.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                data = fd[1:]
                ot.write(fd[0])
                all = 0 
                for i in data:
                    all += float(i)
                for i in data:
                    num = round(float(i)/all,9)
                    ot.write("\t{}".format(num))
                ot.write("\n")
        self.end()
