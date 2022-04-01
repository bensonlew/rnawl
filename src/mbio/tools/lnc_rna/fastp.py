#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob
import re


class FastpAgent(Agent):
    """
    seqprep:用于对SE序列Phred64转Phred33的工具
    author: shicaiping
    """

    def __init__(self, parent):
        super(FastpAgent, self).__init__(parent)
        options = [
            {"name": "in", "type": "infile", "format": "sequence.fastq"},  # 输入文件SE序列
            {"name": "out", "type": "outfile", "format": "sequence.fastq"},  # 输出文件SE序列
        ]
        self.add_option(options)
        self.step.add_steps('fastp')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fastp.start()
        self.step.update()

    def step_end(self):
        self.step.fastp.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("in").is_set:
            raise OptionError("请传入SE序列文件", code = "33705101")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 9
        self._memory = '10G'

    def end(self):
        super(FastpAgent, self).end()


class FastpTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(FastpTool, self).__init__(config)
        self.fastp_path = 'bioinfo/seq/fastp'

    def fastp(self):
        fq_s_path = self.option("in").prop['path']
        fq_s_out = self.output_dir + "/" + os.path.basename(fq_s_path)
        cmd = self.fastp_path + ' -i {} -o {} -A -6 -Q -L -w 8'.format(fq_s_path, fq_s_out)
        self.logger.info(cmd)
        self.logger.info("开始运行fastp")
        command = self.add_command("fastp", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行fastp完成")
        else:
            self.set_error("运行fastp运行出错!", code = "33705102")
            return False

    def run(self):
        """
        运行
        """
        super(FastpTool, self).run()
        self.fastp()
        self.end()
