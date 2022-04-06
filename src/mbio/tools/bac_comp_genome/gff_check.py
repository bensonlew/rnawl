# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.10.14

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file
import os


class GffCheckAgent(Agent):
    """
    主要处理gff文件的处理
    """

    def __init__(self, parent):
        super(GffCheckAgent, self).__init__(parent)
        options = [
            {"name": "gff", "type": "string"},  # gff文件
            {"name": "fa", "type": "infile", "format": "sequence.fasta"}  #gff对应的序列文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gff"):
            raise OptionError("必须输入gff文件！")
        if not self.option("fa").is_set:
            raise OptionError("必须输入fa文件！")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GffCheckAgent, self).end()


class GffCheckTool(Tool):
    def __init__(self, config):
        super(GffCheckTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "/miniconda2/bin/python"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.gff = self.option("gff")
        self.fa = self.option("fa").prop['path']
        prefix = os.path.basename(self.gff).split(".gff")[0]
        self.out = self.output_dir + "/" + str(prefix)

    def run(self):
        """
        运行
        :return:
        """
        super(GffCheckTool, self).run()
        self.run_get_fasta()
        self.end()

    def run_get_fasta(self):
        cmd = '{} {}gff_check.py -g {} -s {} -o {}'.format(self.python_path, self.package_path, self.gff, self.fa, self.out)
        command = self.add_command("run_get_fasta", cmd).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_get_fasta运行完成！")
        else:
            self.set_error("run_get_fasta运行完成运行出错!")
