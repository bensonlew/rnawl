# -*- coding: utf-8 -*-
#__author__ = 'gao.hao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import shutil
import os

class SplitFileAgent(Agent):
    def __init__(self, parent):
        super(SplitFileAgent, self).__init__(parent)
        options = [
            {"name": "files", "type": "infile", "format": "sequence.profile_table"},
            {"name": "lines", "type": "int", "default": 500000},  #行数
            {"name": "out_dir", "type": "outfile", "format": "annotation.mg_anno_dir"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("files").is_set:
            raise OptionError("请传入files文件")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

class SplitFileTool(Tool):
    """
    version 1.0对fastq序列进行拆分
    """
    def __init__(self, config):
        super(SplitFileTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + "/metagbin/"

    def split_fastq(self):
        """
        切分fastq文件
        :return:
        """
        lines = self.option('lines')
        infile = self.option('files').prop['path']
        prefix = os.path.basename(self.option('files').prop['path'])
        sample = self.option('sample')
        if os.path.exists(self.work_dir + "/dir"):
            shutil.rmtree(self.work_dir + "/dir")
        os.mkdir(self.work_dir + "/dir")
        cmd = "{}split_reads.sh {} {} {}".format(self.sh_path, lines, infile, self.work_dir + "/dir/" + prefix)
        command = self.add_command("run_split", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("设置切分结果目录成功!")
            self.set_output()
        else:
            self.set_error("设置切分结果目录运行出错!")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if os.path.exists(self.output_dir + "/dir"):
            shutil.rmtree(self.output_dir + "/dir")
        shutil.copytree(self.work_dir + "/dir",self.output_dir)
        self.option("out_dir", self.output_dir + "/dir")

    def run(self):
        super(SplitFileTool, self).run()
        self.split_fastq()
        self.end()