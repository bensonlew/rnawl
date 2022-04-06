#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bacgenome.common import sum_stat


class BlastDetailAgent(Agent):
    """
    细菌基因组比对结果整理
    version 1.0
    author: gaohao
    last_modify: 2018.04.13
    """

    def __init__(self, parent):
        super(BlastDetailAgent, self).__init__(parent)
        options = [
            {"name": "stat_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "faa_file", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample", "type": "string"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('faa_file').is_set:
            raise OptionError("请设置样本序列文件文件！")
        if not self.option('stat_file').is_set:
            raise OptionError("请设置样本比对结果的结果文件！")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BlastDetailAgent, self).end()


class BlastDetailTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BlastDetailTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/blast_detail.py'

    def run_stat(self):
        sample = self.option('sample')
        faa = self.option('faa_file').prop["path"]
        stat_file = self.option('stat_file').prop["path"]
        cmd = '{} {} -i {} -fa {} -s {} -o {}'.format(self.python_path, self.script,stat_file,faa,sample,self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("run_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_stat运行完成")
        else:
            self.set_error("run_stat运行出错!")


    def set_output(self):
        self.end()

    def run(self):
        """
        运行
        """
        super(BlastDetailTool, self).run()
        self.run_stat()
        self.set_output()


