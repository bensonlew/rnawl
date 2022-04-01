# -*- coding: utf-8 -*-
#__author__ = 'gaohao'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import shutil
import os

class SplitSeqAgent(Agent):
    def __init__(self, parent):
        super(SplitSeqAgent, self).__init__(parent)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "in_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "read_type", "type": "int", "default": 1},
            {"name": "sample", "type": "string"},
            {"name": "lines", "type": "int", "default": 4000000},  # 序列数
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("in_fastq").is_set and not self.option("in_fasta").is_set:
            raise OptionError("请传入in_fastq或in_fasta序列文件！", code="")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="")
        if self.option("in_fastq").is_set:
            if self.option('read_type') not in [1, 2]:
                raise OptionError("read_type类型只能为1或2！", code="")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(SplitSeqAgent, self).end()


class SplitSeqTool(Tool):
    """
    version 1.0对fastq序列进行拆分
    """
    def __init__(self, config):
        super(SplitSeqTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + "/metagbin/"

    def split_fastq(self):
        """
        切分fastq文件
        :return:
        """
        lines = self.option('lines')
        infile = self.option('in_fastq').prop['path']
        sample = self.option('sample')
        if os.path.exists(self.output_dir + "/" + sample + "_fastq_" + str(self.option("read_type"))):
            shutil.rmtree(self.output_dir + "/" + sample + "_fastq_" + str(self.option("read_type")))
        os.mkdir(self.output_dir + "/" + sample + "_fastq_" + str(self.option("read_type")))
        cmd = "{}split_reads.sh {} {} {}".format(self.sh_path, lines, infile, self.output_dir + "/" + sample + "_fastq_" + str(self.option("read_type")) + "/" + sample +"_fastq_")
        command = self.add_command("run_split_fastq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("设置切分结果目录成功!")
        else:
            self.set_error("设置切分结果目录运行出错!")

    def split_fasta(self):
        """

        """
        lines = self.option('lines')
        infile = self.option('in_fasta').prop['path']
        sample = self.option('sample')
        if os.path.exists(self.output_dir + "/" + sample + "_fasta/"):
            shutil.rmtree(self.output_dir + "/" + sample + "_fasta/")
        os.mkdir(self.output_dir + "/" + sample + "_fasta/")
        cmd = "{}split_reads.sh {} {} {}".format(self.sh_path, lines, infile,
                                                 self.output_dir + "/" + sample + "_fasta/" + sample + "_fasta")
        command = self.add_command("run_split_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("设置切分结果目录成功!")
        else:
            self.set_error("设置切分结果目录运行出错!")


    def run(self):
        super(SplitSeqTool, self).run()
        if self.option("in_fastq").is_set:
            self.split_fastq()
        elif self.option("in_fasta").is_set:
            self.split_fasta()
        self.end()
