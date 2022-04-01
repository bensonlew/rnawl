# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.05.03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class MakeblastdbAgent(Agent):
    """
    软件:samtools
    """
    def __init__(self, parent):
        super(MakeblastdbAgent, self).__init__(parent)
        options = [
            {"name": "pop_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "dbtype", "type": "string", 'default': "nucl"},
            {"name": "parse_seqids", "type": "bool", "default": False},  # 是否在makeblastdb的时候增加参数-parse_seqids
            {"name": "db_name", "type": "string", "default": "dbname"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("pop_fa"):
            raise OptionError("必须输入pop.fa文件", code="34503901")
        if self.option("dbtype") not in ["nucl", "prot"]:
            raise OptionError("%s：建库类型不合法！",variables=(self.option("dbtype")), code="34503902")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(MakeblastdbAgent, self).end()


class MakeblastdbTool(Tool):
    def __init__(self, config):
        super(MakeblastdbTool,self).__init__(config)
        self.makeblastdb_path = 'bioinfo/align/ncbi-blast-2.3.0+/bin/makeblastdb'

    def run_makeblastdb(self):
        """
        makeblastdb
        """
        cmd = "{} -in {} -dbtype {} -out {}".format(self.makeblastdb_path, self.option("pop_fa").prop["path"],
                                                    self.option("dbtype"),
                                                    self.output_dir + "/{}".format(self.option("db_name")))
        if self.option("parse_seqids"):
            cmd += " -parse_seqids"
        self.logger.info(cmd)
        command = self.add_command("makeblastdb", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("makeblastdb运行完成")
        else:
            self.set_error("makeblastdb运行失败", code="34503901")

    def run(self):
        super(MakeblastdbTool, self).run()
        self.run_makeblastdb()
        self.end()
