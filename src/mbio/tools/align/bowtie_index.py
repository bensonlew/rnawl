# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/25 9:33

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re


class BowtieIndexAgent(Agent):
    '''
    mapsplice 第2步：调用bowtie为劈开的每个fasta文件做index,mapsplice 只用bowtie1做的index,不过本工具可以做成任何软件都可调用的bowtie index工具
    author: linfang.jin
    last_modify: 2017/1/25 9:33
    '''

    def __init__(self, parent):
        super(BowtieIndexAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", 'format': 'sequence.fasta'},
            {"name": "bowtie_version", "type": "string", "default": "bowtie2"}
        ]
        self.add_option(options)
        self.step.add_steps('bowtie_index')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bowtie_index.start()
        self.step.update()

    def step_end(self):
        self.step.bowtie_index.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not self.option('ref_fa'):
            raise OptionError("必须设置输入的fasta文件")
        if self.option('bowtie_version') not in ('bowtie1', 'bowtie2'):
            raise OptionError("bowtie_version只能被设定为'bowtie1'或'bowtie2'")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        """
        agent结束后一些文件德操作

        :return:"""
        super(BowtieIndexAgent, self).end()


class BowtieIndexTool(Tool):
    def __init__(self, config):
        super(BowtieIndexTool, self).__init__(config)
        self._version = "v2.1.8"
        self.cmd_path_dic = {"bowtie1": "bioinfo/rna/MapSplice-v2.1.8/bin/bowtie-build",
                             "bowtie2": "bioinfo/align/bowtie2-2.2.9/bowtie2-build"} #modify by qingchen.zhang

    def run_bowtie_index(self):
        """
        运行bowtie_index
        :return:
        """

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        cmd = "{} {} {}/{}".format(self.cmd_path_dic[self.option('bowtie_version')],
                                                 self.option('ref_fa').path, self.output_dir,
                                                 os.path.basename(self.option('ref_fa').path))  #modify by qingchen.zhang@20190125
        self.logger.info('运行bowtie_index')
        command = self.add_command("bowtie_index_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bowtie_index运行完成")
        else:
            self.set_error("bowtie_index运行出错!")

    def run(self):
        """
        运行bowtie_index
         :return:
        """
        super(BowtieIndexTool, self).run()
        self.run_bowtie_index()
        self.end()
