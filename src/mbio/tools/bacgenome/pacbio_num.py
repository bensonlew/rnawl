# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2020.06.30

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class PacbioNumAgent(Agent):
    """
    细菌基因组三代数据canu质控后fasta进行统计区间reads数量统计
    """
    def __init__(self, parent):
        super(PacbioNumAgent, self).__init__(parent)
        options = [
            {"name": "input_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"}  # 样品名称
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_fa").is_set:
            raise OptionError("必须添加样品的input_fa的三代fa文件！", code="32101401")
        if not self.option("sample_name"):
            raise OptionError("必须添加样品名称！", code="32101402")

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(PacbioNumAgent, self).end()

class PacbioNumTool(Tool):
    def __init__(self, config):
        super(PacbioNumTool, self).__init__(config)
        self.fq = self.option("input_fa").prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.out_file = self.work_dir + '/pacbio_qc'
        self.pacbio_summary = self.work_dir + '/pacbio.stat.xls'


    def run_num_graphic(self):
        cmd = '{} {}pacbio_num.pl {} 1000 1000 {}'.format(self.perl_path,self.perl_script,self.fq,self.option("sample_name"))
        command = self.add_command("run_num_graphic", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_num_graphic运行完成")
        else:
            self.set_error("run_num_graphic运行出错!", code="32101401")

    def set_output(self):
        for file in [self.option("sample_name") + '.clean.len.xls',self.option("sample_name") + '.PacBio_statistics.xls']:
            if os.path.exists(self.output_dir + '/' + file):
                os.remove(self.output_dir + '/' + file)
            os.link(self.work_dir + '/' + file,self.output_dir + '/' + file)

    def run(self):
        super(PacbioNumTool, self).run()
        self.run_num_graphic()
        self.set_output()
        self.end()