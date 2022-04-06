#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class MobileStatAgent(Agent):
    """
    细菌可移动元件统计
    version 1.0
    author: gaohao
    last_modify: 2020.10.14
    """
    def __init__(self, parent):
        super(MobileStatAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('dir'):
            raise OptionError("请设置dir文件目录！")
        if not self.option('genome').is_set:
            raise OptionError("请设置gemone序列文件！")
        if not self.option('sample_name'):
            raise OptionError("请设置样品名！")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(MobileStatAgent, self).end()


class MobileStatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MobileStatTool, self).__init__(config)
        self.scaf_seq = self.option("genome").prop["path"]
        self.dir =self.option('dir')
        self.sample = self.option('sample_name')
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + "/mobile_genetic_elements/"

    def run_mobilestat(self):
        cmd = '{} {}mobile_stat.py --g {} --n {} --d {} --o {}'.format(self.python_path, self.script, self.scaf_seq, self.sample, self.dir, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("run_mobilestat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_mobilestat运行完成")
        else:
            self.set_error("run_mobilestat运行出错!")

    def run(self):
        """
        运行
        """
        super(MobileStatTool, self).run()
        self.run_mobilestat()
        self.end()
