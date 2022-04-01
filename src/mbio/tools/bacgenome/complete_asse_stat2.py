# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/9'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file,link_dir


class CompleteAsseStat2Agent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(CompleteAsseStat2Agent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta", "required": True},
            # {"name": "genome_id", "type": "string", "required": True},
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"


class CompleteAsseStat2Tool(Tool):
    def __init__(self, config):
        super(CompleteAsseStat2Tool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.python_path = "program/Python/bin/python"
        # self.fa = os.path.join(self.work_dir, "all.fa")
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/complete_get_seq.py"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/complete_stat2.pl"

    def run_asse_seq(self):
        """
        description
        :return:
        """
        # cmd = '{} {} {}'.format(self.perl_path, self.perl_script, self.fa)
        cmd = '{} {} {} {}'.format(self.perl_path, self.perl_script, self.option("fa").prop["path"], self.option("sample_name"))
        command = self.add_command("asse_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("asse_seq运行完成")
        else:
            self.set_error("asse_seq运行出错！")

    def set_output(self):
        # self.option("out_fa").set_path(self.fa)
        stat1_name = "assemble.stat.xls"
        stat2_name = "assemble.summary.xls"
        link_file(os.path.join(self.work_dir, stat1_name), os.path.join(self.output_dir, stat1_name))
        link_file(os.path.join(self.work_dir, stat2_name), os.path.join(self.output_dir, stat2_name))
        # link_file(self.fa, os.path.join(self.output_dir, "all.fa"))
        link_dir(os.path.join(self.work_dir, "seq_dir"), os.path.join(self.output_dir, "seq_dir"))

    def run(self):
        super(CompleteAsseStat2Tool, self).run()
        # self.get_seq()
        self.run_asse_seq()
        self.set_output()
        self.end()