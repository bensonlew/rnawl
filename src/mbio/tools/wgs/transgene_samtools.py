# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.04.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class TransgeneSamtoolsAgent(Agent):
    """
    软件:samtools
    """
    def __init__(self, parent):
        super(TransgeneSamtoolsAgent, self).__init__(parent)
        options = [
            {"name": "pop_fa", "type": "infile", "format": "sequence.fasta"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pop_fa"):
            raise OptionError("必须输入pop.fa文件", code="34507101")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(TransgeneSamtoolsAgent,self).end()


class TransgeneSamtoolsTool(Tool):
    def __init__(self, config):
        super(TransgeneSamtoolsTool, self).__init__(config)
        self.samtools_path = 'miniconda2/bin/samtools'
        self.makeblastdb_path = 'bioinfo/align/ncbi-blast-2.3.0+/bin/makeblastdb'

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
            self.set_error("samtools运行失败", code="34507101")

    def run_makeblastdb(self):
        """
        makeblastdb
        """
        cmd = "{} -in {} -dbtype nucl -out {}".format(self.makeblastdb_path, self.option("pop_fa").prop["path"],
                                                      self.output_dir + "/dbname")
        self.logger.info(cmd)
        command = self.add_command("makeblastdb", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("makeblastdb运行完成")
        else:
            self.set_error("makeblastdb运行失败", code="34507102")

    def run(self):
        super(TransgeneSamtoolsTool, self).run()
        self.run_faidx()
        self.run_makeblastdb()
        self.end()
