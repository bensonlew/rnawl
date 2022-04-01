# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import pandas as pd
import subprocess, os,re


class PlotTreeAgent(Agent):
    """
    比较基因组16s进化树分析
    """

    def __init__(self, parent):
        super(PlotTreeAgent, self).__init__(parent)
        options = [
            {"name": "16s_dir", "type": "infile", "formate": "sequence.fasta_dir"},  # 16s的
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("16s_dir"):
            raise OptionError("必须输入16s_dir文件夹！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(PlotTreeAgent, self).end()


class PlotTreeTool(Tool):
    def __init__(self, config):
        super(PlotTreeTool, self).__init__(config)
        self._version = "1.0"
        self.python_path ="/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/common_anno.py'
        self.databases = self.option("database_list").split(",")
        self.dir = self.option("dir")

    def run(self):
        """
        运行
        :return:
        """
        super(PlotTreeTool, self).run()
        self.run_anno()
        self.end()

    def run_anno(self):
        for i in self.databases:
            if not os.path.exists(self.output_dir + "/" + i):
                os.mkdir(self.output_dir + "/" + i)
            cmd = '{} {} -dir {} -t {} -o {}'.format(self.python_path, self.python_script, self.dir, i, self.output_dir + "/" + i + "/")
            command = self.add_command("run_{}".format(i), cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_{}运行完成！".format(i))
            else:
                self.set_error("run_{}运行完成运行出错!".format(i))