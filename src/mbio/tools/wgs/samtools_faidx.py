# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.05.03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os

class SamtoolsFaidxAgent(Agent):
    """
    软件:samtools
    """
    def __init__(self, parent):
        super(SamtoolsFaidxAgent, self).__init__(parent)
        options = [
            {"name": "pop_fa", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pop_fa"):
            raise OptionError("必须输入pop.fa文件", code="34505001")

    def set_resource(self):
        self._cpu = 2
        self._memory = "3G"

    def end(self):
        super(SamtoolsFaidxAgent,self).end()

class SamtoolsFaidxTool(Tool):
    def __init__(self, config):
        super(SamtoolsFaidxTool, self).__init__(config)
        self.samtools_path = 'miniconda2/bin/samtools'

    def run_faidx(self):
        """
        samtools faidx
        """
        cmd = "{} faidx {}".format(self.samtools_path, self.option("pop_fa").prop["path"])
        self.logger.info(cmd)
        command = self.add_command("samtools_faidx", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools运行完成")
        else:
            self.set_error("samtools运行失败", code="34505001")

    def run(self):
        super(SamtoolsFaidxTool, self).run()
        self.run_faidx()
        self.end()