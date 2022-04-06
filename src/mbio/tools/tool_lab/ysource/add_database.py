# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import datetime
import random
import re
import os
import sys
import shutil

class AddDatabaseAgent(Agent):
    """
    把分型的结果导入到 Mysql数据库里
    """
    def __init__(self, parent):
        super(AddDatabaseAgent, self).__init__(parent)
        options = [ 
            {"name":"source_result","type":"string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        if not os.path.exists(self.option("source_result")):
            raise OptionError("请检查分型结果是否生成")
    
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AddDatabaseAgent, self).end()
    
class AddDatabaseTool(Tool):

    def __init__(self, config):
        super(AddDatabaseTool, self).__init__(config)
        self._version = '1.0'
        self.python = 'miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + '/tool_lab/yoogene/yfull/add_source_database.py'

    def run(self):
        super(AddDatabaseTool, self).run()
        self.run_adddatabase()
        self.end()
        
    def run_adddatabase(self):
        """
        python add_sourse_database.py -i inputfile
        """
        cmd = "{} {} ".format(self.python,self.script)
        cmd += "-i {}".format(self.option("source_result"))
        self.logger.info(cmd)
        self.logger.info("开始计算")
        command1 = self.add_command("add_database",cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("导入MySQL运行完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
