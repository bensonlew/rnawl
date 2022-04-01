# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2021.04.08

import os, time
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir
import json

class MlstTreeAgent(Agent):
    """
    多个基因组MLST进化树分析
    """
    def __init__(self, parent):
        super(MlstTreeAgent, self).__init__(parent)
        options = [
            {"name": "martix", "type": "infile", "format": "sequence.profile_table"},# ST型文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("martix").is_set:
            raise OptionError("必须输入martix文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(MlstTreeAgent, self).end()

class MlstTreeTool(Tool):
    def __init__(self, config):
        super(MlstTreeTool, self).__init__(config)
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "tool_lab/mlst_tree.sh")
        self.script = self.config.PACKAGE_DIR + '/tool_lab/martix_calculate.py'
        self.path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/rapidNJ-master/bin"
        self.set_environ(PATH=self.path)

    def run_mlsttree(self):
        """
        运行mlst
        :return:
        """
        cmd = "{} {} {} {} ".format(self.shell, self.shell_path, self.script, self.option("martix").prop['path'])
        self.logger.info(cmd)
        command = self.add_command("run_mlsttree", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_mlsttree运行完成！")
        else:
            self.set_error("run_mlsttree运行出错!")

    def set_output(self):
        """
        设置结果文件目录
        """
        os.link(self.work_dir+"/mlst_tree.nwk", self.output_dir+"/mlst_tree.nwk")
        self.logger.info("生成结果文件完成")

    def run(self):
        """
        运行
        """
        super(MlstTreeTool, self).run()
        self.run_mlsttree()
        self.set_output()
        self.end()