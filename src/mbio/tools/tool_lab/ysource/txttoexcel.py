# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from mbio.packages.whole_transcriptome.utils import runcmd
# import pandas as pd
import unittest
import xlsxwriter
import datetime
import random
import re
import os
import sys
import shutil

class TxttoexcelAgent(Agent):
    def __init__(self, parent):
        super(TxttoexcelAgent, self).__init__(parent)
        self.options = [
            {"name":"input_file","type":"string"},
            {"name":"output_name","type":"string"},
        ]
        self.add_option(self.options)

    def check_option(self):
        """
        参数检查
        """
        if not os.path.exists(self.option("input_file")):
            raise OptionError("找不到{}".format(option("input_file")))
        if not self.option("output_name"):
            raise OptionError("请设置输出文件名")
    
    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"
    
    def end(self):
        super(TxttoexcelAgent, self).end()

class TxttoexcelTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(TxttoexcelTool, self).__init__(config)
        self._version = "v1.0"
    
    def run(self):
        super(TxttoexcelTool, self).run()
        self.run_txtToExcel(self.option("input_file"),self.option("output_name"))
        self.set_output()
        self.end()

    def run_txtToExcel(self,input, output):
        if os.path.exists(input):
            f = open(input)
            wo = xlsxwriter.Workbook(output)
            sheet = wo.add_worksheet()
            x = 0
            while 1:
                line = f.readline()
                if not line:
                    break
                for i in range(len(line.split('\t'))):
                    item = line.split('\t')[i]
                    sheet.write(x, i, item)
                x += 1
            wo.close()
            f.close()

    def set_output(self):
        """
        设置输出文件夹
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "{}".format(self.option("output_name"))),
                    os.path.join(self.output_dir, "{}".format(self.option("output_name"))))        
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))