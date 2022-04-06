#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os

class CoverageSumAgent(Agent):
    """
    用于计算scaffolds的测序深度的汇总统计
    version 1.0
    author: gaohao
    last_modify: 2019.08.07
    """

    def __init__(self, parent):
        super(CoverageSumAgent, self).__init__(parent)
        options = [
            {"name": "depth_files", "type": "infile", "format": "metagbin.depth_dir"},
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('depth_files'):
            raise OptionError("reads mapping的depth_files文件字段不存在！")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(CoverageSumAgent, self).end()


class CoverageSumTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(CoverageSumTool, self).__init__(config)
        self.depth =self.option('depth_files').prop['path']
        self.python_path = "miniconda2/bin/python"
        self.cover =self.config.PACKAGE_DIR + "/metagbin/metagbin_coverage.py"

    def run_summarize_depth(self):
        cmd = "{} {} {} {}".format(self.python_path,self.cover, self.depth, self.work_dir + "/all.depth.txt")
        self.logger.info(cmd)
        self.logger.info("开始运行metagbin_coverage")
        command = self.add_command("metagbin_coverage", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行metagbin_coverage完成")
        else:
            self.set_error("运行metagbin_coverage运行出错!")

    def set_output(self):
        if os.path.exists(self.output_dir + "/depth.txt"):
            os.remove(self.output_dir + "/depth.txt")
        os.link(self.work_dir + "/all.depth.txt", self.output_dir + "/depth.txt")

    def run(self):
        """
        运行
        """
        super(CoverageSumTool, self).run()
        self.run_summarize_depth()
        self.set_output()
        self.end()