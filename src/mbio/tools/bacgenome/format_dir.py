#-*- coding: utf-8 -*
import os
import re, shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class FormatDirAgent(Agent):
    """
    处理format文件，处理成Integron_finder程序所需要的输入文件
    目的是根据基因组文件、gff文件和序列文件进行format
    """
    def __init__(self, parent):
        super(FormatDirAgent, self).__init__(parent)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  # 细菌的基因组fasta文件或宏基因组组装的fasta文件
            {"name": "gene_faa", "type": "infile",  "format": "sequence.fasta"},  # 预测的基因核酸文件
            {"name": "gene_gff", "type": "string"},   # 预测的基因gff文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome_fa").is_set:
            raise OptionError("必须设置参数genome_fa文件!")
        if not self.option("gene_faa").is_set:
            raise OptionError("必须设置参数gene_faa文件!")
        if not self.option("gene_gff"):
            raise OptionError("必须设置参数gene_gff文件!")

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(FormatDirAgent, self).end()

class FormatDirTool(Tool):
    def __init__(self, config):
        super(FormatDirTool, self).__init__(config)
        self.python = "/program/Python/bin/python"
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/format_dir.py"
        self.genome = self.option("genome_fa").prop['path']
        self.gene_faa =self.option("gene_faa").prop['path']
        self.gene_gff =self.option("gene_gff")

    def run_format(self):
        """
        运行脚本进行format文件
        :return:
        """
        cmd = "{} {} --s {} --f {} --g {} --o {}".format(self.python, self.python_script, self.genome, self.gene_faa, self.gene_gff, self.output_dir)
        command = self.add_command("run_format", cmd).run()
        self.logger.info(cmd)
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_format运行完成！")
        else:
            self.set_error("run_format运行完成运行出错!")

    def run(self):
        super(FormatDirTool, self).run()
        self.run_format()
        self.end()
