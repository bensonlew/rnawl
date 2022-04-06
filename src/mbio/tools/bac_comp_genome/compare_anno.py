# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.09.27

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
import pandas as pd
import subprocess, os,re


class CompareAnnoAgent(Agent):
    """
    处理多个样品的注释汇总
    """

    def __init__(self, parent):
        super(CompareAnnoAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},  # 结果目录文件夹
            {"name": "database_list", "type": "string", "default": "cog,kegg,cazy,card,vfdb,tcdb,phi,secretory,tmhmm,signalp,all"},  # 基因具体的注释信息
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("dir"):
            raise OptionError("必须注释输入文件！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(CompareAnnoAgent, self).end()


class CompareAnnoTool(Tool):
    def __init__(self, config):
        super(CompareAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.python_path ="/miniconda2/bin/python"
        self.python_script = self.config.PACKAGE_DIR + '/bac_comp_genome/common_anno.py'
        self.databases = self.option("database_list").split(",")
        self.dir = self.option("dir")

    def run(self):
        """
        运行
        :return:
        """
        super(CompareAnnoTool, self).run()
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