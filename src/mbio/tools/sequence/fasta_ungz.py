# -*- coding: utf-8 -*-
# __author__ :gaohao
# last_modify: 20190116

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os,re
import math


class FastaUngzAgent(Agent):
    def __init__(self, parent):
        super(FastaUngzAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "string"},
            {"name": "sample_name", "type": "string"},
            {"name": "out_fa", "type": "outfile","format":"sequence.fasta"},
        ]
        self.add_option(options)
        self.step.add_steps("fasta_ungz")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fasta_ungz.start()
        self.step.update()

    def stepfinish(self):
        self.step.fasta_ungz.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fasta'):
            raise OptionError("必须输入fasta文件")
        if not self.option('sample_name'):
            raise OptionError("必须输入样品名")
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        memory = os.path.getsize(self.option('fasta').split(" ")[0])
        n = memory / (1024 * 1024 * 1024)
        n = math.ceil(n)
        self._memory = '{}G'.format(int(n))

    def end(self):
        super(FastaUngzAgent, self).end()


class FastaUngzTool(Tool):
    def __init__(self, config):
        super(FastaUngzTool, self).__init__(config)
        self._version = "v1.0"
        self.sh_path =  self.config.PACKAGE_DIR +"/sequence/scripts/unzip.sh"

    def run(self):
        super(FastaUngzTool, self).run()
        self.run_unzip()
        self.set_output()
        self.end()

    def run_unzip(self):
        sample = self.option("sample_name")
        self.fasta = self.option("fasta")
        f_name = sample + ".assemble"  + ".fa"
        new_fasta = os.path.join(self.work_dir, f_name)
        if re.search(r'.tar.gz', self.fasta):
            unzip_cmd = "tar -zxvf " + self.fasta
        elif re.search(r'.gz', self.fasta):
            unzip_cmd = "zcat " + self.fasta + " > " + new_fasta
        self.logger.info("start unzip")
        self.logger.info("command: %s" % unzip_cmd)
        try:
            subprocess.check_output(unzip_cmd, shell=True)
            self.logger.info("unzip done")
        except subprocess.CalledProcessError:
            self.set_error("unzip error")
            raise Exception("unzip error")

    def set_output(self):
        if re.search(r'.tar.gz', self.fasta):
            n = re.search("(.*).tar.gz", self.fasta.split("/")[-1])
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + ".assemble.fa"):
                os.remove(self.output_dir + '/' + self.option("sample_name") + ".assemble.fa")
            os.link(self.work_dir + '/' + n.group(1),self.output_dir + '/' + self.option("sample_name") + ".assemble.fa")
            self.option('out_fa', self.output_dir + '/' + self.option("sample_name") + ".assemble.fa")
        elif re.search(r'.gz', self.fasta):
            if os.path.exists(self.output_dir + '/' + self.option("sample_name") + ".assemble.fa"):
                os.remove(self.output_dir + '/' + self.option("sample_name") + ".assemble.fa")
            os.link(self.work_dir + '/' + self.option("sample_name") + ".assemble.fa",
                    self.output_dir + '/' + self.option("sample_name") + ".assemble.fa")
            self.option('out_fa', self.output_dir + '/' + self.option("sample_name") + ".assemble.fa")