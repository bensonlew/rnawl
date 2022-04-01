# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/5/9'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import link_file,link_dir


class CompleteGetSeqAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(CompleteGetSeqAgent, self).__init__(parent)
        options = [
            {"name": "fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "genome_id", "type": "string", "required": True},
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


class CompleteGetSeqTool(Tool):
    def __init__(self, config):
        super(CompleteGetSeqTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.fa = os.path.join(self.work_dir, "all.fa")
        self.python_script = self.config.PACKAGE_DIR + "/bacgenome/complete_get_seq.py"

    def get_seq(self):
        cmd = "{} {} --f {} --id {} --o {}".format(self.python_path, self.python_script, self.option("fa").prop["path"], self.option("genome_id"), self.fa)
        if self.config.DBVersion:
            cmd += ' --m ' + str(self.config.DBVersion)
        command = self.add_command("get_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("seq运行完成")
        else:
            self.set_error("seq运行出错！")

    def set_output(self):
        self.option("out_fa").set_path(self.fa)

    def run(self):
        super(CompleteGetSeqTool, self).run()
        self.get_seq()
        self.set_output()
        self.end()