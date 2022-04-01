# -*- coding: utf-8 -*-
# __author__ :zhaobinbin
# last_modify: 20200618

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.metagbin.common_function import link_dir
import os
import math
import shutil


class FastqUngzAgent(Agent):
    def __init__(self, parent):
        super(FastqUngzAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},
        ]
        self.add_option(options)


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError("必须输入fastq.gz文件")
        if not self.option('fastq2'):
            raise OptionError("必须输入fastq.gz文件")
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        memory = os.path.getsize(self.option('fastq1').prop["path"].split(" ")[0])*2
        n = memory / (1024 * 1024 * 1024)
        n = math.ceil(n)
        self._memory = '{}G'.format(int(n))

    def end(self):
        super(FastqUngzAgent, self).end()


class FastqUngzTool(Tool):
    def __init__(self, config):
        super(FastqUngzTool, self).__init__(config)
        self._version = "v1.0"
        self.sh_path = self.config.PACKAGE_DIR +"/sequence/scripts/unzip.sh"

    def run(self):
        super(FastqUngzTool, self).run()
        self.run_unzip(1)
        self.run_unzip(2)
        self.set_output()
        self.end()

    def run_unzip(self, index):
        if os.path.exists(self.work_dir + "/unzip_dir"):
            pass
        else:
            os.mkdir(self.work_dir + "/unzip_dir")
        path1 = self.option('fastq1').prop["path"]
        path2 = self.option('fastq2').prop["path"]
        self.path1 = os.path.dirname(path1)
        self.path2 = os.path.dirname(path2)
        self.file1 = os.path.basename(path1)
        self.file2 = os.path.basename(path2)
        file_path1 = os.path.join(self.path1, self.file1)
        file_path2 = os.path.join(self.path2, self.file2)
        file_name1 = self.work_dir + "/unzip_dir/" + ".".join(self.file1.split('.')[0:-2])
        file_name2 = self.work_dir + "/unzip_dir/" + ".".join(self.file2.split('.')[0:-2])
        unzip_cmd1 = "zcat " +  file_path1 + " > " + file_name1 + ".fq"
        unzip_cmd2 = "zcat " +  file_path2 + " > " + file_name2 + ".fq"

        self.logger.info("start unzip")
        self.logger.info("command: %s" % unzip_cmd1)
        self.logger.info("command: %s" % unzip_cmd2)
        if index == 1:
            try:
                subprocess.check_output(unzip_cmd1, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error")
                raise Exception("unzip error")
        elif index == 2:
            try:
                subprocess.check_output(unzip_cmd2, shell=True)
                self.logger.info("unzip done")
            except subprocess.CalledProcessError:
                self.set_error("unzip error")
                raise Exception("unzip error")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        link_dir(self.work_dir + "/unzip_dir", self.output_dir)
