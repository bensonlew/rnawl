# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.18

import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class TransposonCombinAgent(Agent):
    """
    主要IS和复杂转座子，兼容处理
    """

    def __init__(self, parent):
        super(TransposonCombinAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},  # 结果文件目录
            {"name": "sample", "type": "string"},  #样品名称
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("dir"):
            raise OptionError("必须设置参数dir路径!")

    def set_resource(self):
        self._cpu = 4
        self._memory = '10G'

    def end(self):
        super(TransposonCombinAgent, self).end()


class TransposonCombinTool(Tool):
    def __init__(self, config):
        super(TransposonCombinTool, self).__init__(config)
        self.python_path = '/miniconda2/bin/python'
        self.script_path = self.config.PACKAGE_DIR + '/mobile_genetic_elements/'
        self.sample = self.option("sample")

    def run_transposon(self):
        files = os.listdir(self.option("dir"))
        if len(files) >=2:
            for file in files:
                if re.search('.is.xls', file):
                    is_path = self.option("dir") + "/" + file
                elif re.search('.transposon.xls', file):
                    trans_path = self.option("dir") + "/" + file
            cmd2 = "{} {}combin_transposon.py --p {} --i {} --o {} ".format(self.python_path, self.script_path, trans_path, is_path, self.output_dir + "/" + self.sample + ".transposon.xls")
            command2 = self.add_command('run_transposon', cmd2).run()
            self.wait(command2)
            self.logger.info(command2)
            if command2.return_code == 0:
                self.logger.info("run_transposon处理成功")
            else:
                self.set_error("run_transposon文件处理失败")
        elif len(files) == 1:
            for file in files:
                os.link(self.option("dir") + "/" + file, self.output_dir + "/" +self.sample + ".transposon.xls")


    def run(self):
        super(TransposonCombinTool, self).run()
        self.run_transposon()
        self.end()
