# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0521

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BwaConfigAgent(Agent):
    """
    软件:bwa index
    """
    def __init__(self, parent):
        super(BwaConfigAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("必须输入ref.fa文件", code="34500901")

    def set_resource(self):
        self._cpu = 2
        self._memory = "50G"

    def end(self):
        super(BwaConfigAgent, self).end()


class BwaConfigTool(Tool):
    def __init__(self, config):
        super(BwaConfigTool, self).__init__(config)
        self.bwa_path = 'bioinfo/align/bwa-0.7.15/bwa'

    def run_bwa(self):
        """
        bwa index
        """
        cmd = "{} index {}".format(self.bwa_path, self.option("ref_fa").prop['path'])
        command = self.add_command("fa_bwa", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bwa运行完成")
        else:
            self.set_error("bwa运行失败", code="34500901")

    def run(self):
        super(BwaConfigTool, self).run()
        self.run_bwa()
        self.end()
