# -*- coding: utf-8 -*-
#__author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import shutil
import os
from mbio.packages.metagbin.common_function import link_dir
from mbio.packages.metagbin.common_function import link_file
class SplitFastqAgent(Agent):
    def __init__(self, parent):
        super(SplitFastqAgent, self).__init__(parent)
        options = [
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq"},
            {"name": "read_type", "type": "string"},
            {"name": "sample", "type": "string"},
            {"name": "lines", "type": "int", "default": 40000000},  # 序列数
        ]
        self.add_option(options)
    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("in_fastq").is_set:
            raise OptionError("请传入in_fastq序列文件", code="")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="")
        if self.option('read_type') not in ['1','2']:
            raise OptionError("read_type类型只能为l或r！", code="")
    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'
class SplitFastqTool(Tool):
    """
    version 1.0对fastq序列进行拆分
    """
    def __init__(self, config):
        super(SplitFastqTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + "/metagbin/"
    def split_fastq(self):
        """
        切分fastq文件
        :return:
        """
        lines = self.option('lines')
        infile = self.option('in_fastq').prop['path']
        sample = self.option('sample')
        if os.path.exists(self.work_dir + "/" + sample + "_fastq_" + self.option("read_type")):
            shutil.rmtree(self.work_dir + "/" + sample + "_fastq_" + self.option("read_type"))
        os.mkdir(self.work_dir + "/" + sample + "_fastq_" + self.option("read_type"))
        cmd = "{}split_reads.sh {} {} {}".format(self.sh_path, lines, infile, self.work_dir + "/" + sample + "_fastq_" + self.option("read_type") + "/" + sample +"_fastq_")
        command = self.add_command("run_split_fastq", cmd).run()
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
        sample = self.option('sample')
        if os.path.exists(self.output_dir + "/" + sample + "_fastq_"+ self.option("read_type")):
            shutil.rmtree(self.output_dir + "/" + sample + "_fastq_"+ self.option("read_type"))
        link_dir(self.work_dir +  "/" + sample + "_fastq_"+  self.option("read_type"),self.output_dir +  "/" + sample + "_fastq_"+ self.option("read_type"))
    def run(self):
        super(SplitFastqTool, self).run()
        self.split_fastq()
        self.end()