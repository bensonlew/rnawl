# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file
import os


class MergeFastaAgent(Agent):
    """
    16s序列进行比对对齐和掩盖低置信度的位置处理
    """

    def __init__(self, parent):
        super(MergeFastaAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},  # fasta文件夹
            {"name": "sample_list", "type": "string"},  # 文件中有样品信息
            {"name": "out", "type": "outfile", "format": "sequence.fasta"}, #比对对齐的序列文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta_dir").is_set:
            raise OptionError("必须输入fasta_dir文件夹！")
        if not self.option("sample_list"):
            raise OptionError("必须输入sample_list！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(MergeFastaAgent, self).end()


class MergeFastaTool(Tool):
    def __init__(self, config):
        super(MergeFastaTool, self).__init__(config)
        self._version = "1.0"
        self.fasta_dir = self.option("fasta_dir").prop["path"]
        self.sample_list = self.option("sample_list")
        self.python_path = "/program/Python/bin/python"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.out = self.output_dir + "/all.align.fasta"

    def run(self):
        """
        运行
        :return:
        """
        super(MergeFastaTool, self).run()
        self.run_merge_fasta()
        self.set_output()
        self.end()

    def run_merge_fasta(self):
        cmd = '{} {}merge_fasta.py -t {} -s {} -o {}'.format(self.python_path, self.package_path, self.fasta_dir, self.sample_list, self.out)
        command = self.add_command("run_merge_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_merge_fasta运行完成！")
        else:
            self.set_error("run_merge_fasta运行完成运行出错!")

    def set_output(self):
        self.option("out", self.output_dir + "/all.align.fasta")