# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.06.08

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class PacbioCleanAgent(Agent):
    """
    细菌基因组三代数据canu质控后fasta进行统计
    """
    def __init__(self, parent):
        super(PacbioCleanAgent, self).__init__(parent)
        options = [
            {"name": "input_fq", "type": "infile", "format": "sequence.fastq"},
            {"name": "type", "type": "string"}, # pacbio or nanopore
            {"name": "sample_name", "type": "string"}  # 样品名称
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_fq").is_set:
            raise OptionError("必须添加样品的input_fq的三代fa文件！", code="32101401")
        if not self.option("sample_name"):
            raise OptionError("必须添加样品名称！", code="32101402")

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(PacbioCleanAgent, self).end()

class PacbioCleanTool(Tool):
    def __init__(self, config):
        super(PacbioCleanTool, self).__init__(config)
        self.fq = self.option("input_fq").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/sequence/"
        self.out_file = self.work_dir + '/pacbio_qc'
        self.pacbio_summary = self.work_dir + '/pacbio.stat.xls'


    def run_len_graphic(self):
        cmd = '{} {}pacbio_length.pl {} 100 1000 {} {}'.format(self.perl_path,self.perl_script,self.fq,self.option("type"), self.option("sample_name"))
        command = self.add_command("run_len_graphic", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_len_graphic运行完成")
        else:
            self.set_error("run_len_graphic运行出错!", code="32101401")

    def set_output(self):
        if self.option("type") == "pacbio":
            for file in [self.option("sample_name") + '.clean.len.xls',
                         self.option("sample_name") + '.PacBio_statistics.xls']:
                if os.path.exists(self.output_dir + '/' + file):
                    os.remove(self.output_dir + '/' + file)
                os.link(self.work_dir + '/' + file, self.output_dir + '/' + file)
        elif self.option("type") == "nanopore":
            for file in [self.option("sample_name") + '.clean.len.xls',
                         self.option("sample_name") + '.Nanopore_statistics.xls']:
                if os.path.exists(self.output_dir + '/' + file):
                    os.remove(self.output_dir + '/' + file)
                os.link(self.work_dir + '/' + file, self.output_dir + '/' + file)


    def run(self):
        super(PacbioCleanTool, self).run()
        self.run_len_graphic()
        self.set_output()
        self.end()