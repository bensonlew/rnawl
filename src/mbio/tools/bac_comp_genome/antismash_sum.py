#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file

class AntismashSumAgent(Agent):
    """
    生成基因组岛文件
    version 1.0
    author: gaohao
    last_modify: 2019.10.11
    """

    def __init__(self, parent):
        super(AntismashSumAgent, self).__init__(parent)
        options = [
            {"name": "dir", "type": "string"},  #
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('dir'):
            raise OptionError("请设置antismash的结果文件夹路径！")


    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(AntismashSumAgent, self).end()


class AntismashSumTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(AntismashSumTool, self).__init__(config)
        self.dir = self.option("dir")
        self.sample_list = self.option("sample_list").prop['path']
        self.python_path = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"

    def run_antismash_stat(self):
        cmd = '{} {}antismash_summary.py -i {} -o {} -s {}'.format(self.python_path, self.python_script, self.dir, self.work_dir, self.sample_list)
        self.logger.info(cmd)
        command = self.add_command("run_antismash_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_antismash_stat运行完成")
        else:
            self.set_error("run_antismash_stat运行出错!")

    def set_output(self):
        link_file(self.work_dir + '/all.antismash_abund.xls', self.output_dir + '/all.antismash_abund.xls')
        link_file(self.work_dir + '/all.antismash_genelist.xls', self.output_dir + '/all.antismash_genelist.xls')

    def run(self):
        """
        运行
        """
        super(AntismashSumTool, self).run()
        self.run_antismash_stat()
        self.set_output()
        self.end()