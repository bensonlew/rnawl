# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/13 14:03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.ref_rna.gene_structure.hdrs import HdrsFile
import subprocess
import re
import shutil
import json
import os
import glob
from mbio.files.sequence.fasta import FastaFile


class AsprofileHdrsAgent(Agent):
    '''
    运行asprofile的count fasta 功能，生成参考基因组序列文件的hdrs文件
    '''

    def __init__(self, parent):
        super(AsprofileHdrsAgent, self).__init__(parent)
        options = [{"name": "ref_fa", "type": "infile", "format": "sequence.fasta", "default": None},
                   ]
        self.add_option(options)
        self.step.add_steps('asprofile_count_fasta')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def check_options(self):
        if self.option('ref_fa') is None:
            raise OptionError("必须设置参考基因组文件")

    def step_start(self):
        self.step.asprofile_count_fasta.start()
        self.step.update()

    def step_end(self):
        self.step.asprofile_count_fasta.finish()
        self.step.update()

    def set_resource(self):
        '''
        所需资源
        :return:
        '''
        self._cpu = 10
        self._memory = '100G'

    def end(self):

        super(AsprofileHdrsAgent, self).end()


class AsprofileHdrsTool(Tool):
    def __init__(self, config):
        super(AsprofileHdrsTool, self).__init__(config)
        self._version = "vb-1.0.4"
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/ASprofile.b-1.0.4/"
        self.python_path = 'program/Python/bin/python'

    def run_hdrs(self):
        """
        :return:
        """
        hdrs_path = os.path.join(self.output_dir,"ref.fa.hdrs")
        cmd = "{}  {}write_hdrs.py  -fa {}  -hdrs  {}".format(self.python_path, self.script_path,
                                                              self.option('ref_fa').path, hdrs_path)
        self.logger.info('开始运行asprofile:count_fasta')
        command = self.add_command("run_hdrs_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("count_fasta运行完成")
        else:
            self.set_error("count_fasta运行出错!")


    def run(self):
        """
        运行asprofile_count_ref，输入文件为fasta格式
         :return:
        """
        super(AsprofileHdrsTool, self).run()
        self.run_hdrs()
        self.end()
