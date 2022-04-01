# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# modified 2018.05.03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BlastNAgent(Agent):
    """
    软件:Blastn
    """

    def __init__(self, parent):
        super(BlastNAgent, self).__init__(parent)
        options = [
            {"name": "query_fa", "type": "infile", "format": "sequence.fasta"},  # 要比对的序列
            {"name": "dbname_nsq", "type": "string"},  # 比对库的路径
            {"name": "outfmt", "type": "string"},  # 输出格式 tab格式是6
            {"name": "num_threads", "type": "int", "default": 8},
            {"name": "evalue", "type": "float", "default": 1e-5},  # 增加参数evalue、num_alignments modified by zengjing 20180508
            {"name": "num_alignments", "type": "int"},
            {"name": "sample_id", "type": "string"}  # 用于标记输出文件名字，可以不输入
        ]
        self.add_option(options)
        self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def check_options(self):
        if not self.option("query_fa"):
            raise OptionError("必须输入query_fa文件", code="34500601")
        if not self.option("dbname_nsq"):
            raise OptionError("必须输入dbname_nsq文件", code="34500602")

    def set_resource(self):
        self._cpu = self.option("num_threads")
        self._memory = "20G"

    def end(self):
        super(BlastNAgent, self).end()


class BlastNTool(Tool):
    def __init__(self, config):
        super(BlastNTool, self).__init__(config)
        self.blastn_path = 'bioinfo/align/ncbi-blast-2.3.0+/bin/blastn'

    def run_blastn(self):
        """Is
        Blastn
        """
        if not self.option("sample_id"):
            outfile_name = self.output_dir + "/blastn.blast"
        else:
            outfile_name = self.output_dir + "/{}.blast".format(self.option("sample_id"))
        if not self.option("outfmt"):
            cmd = "{} -query {} -out {} -db {} -evalue {} -num_threads {}"\
                .format(self.blastn_path, self.option("query_fa").prop["path"],
                        outfile_name, self.option("dbname_nsq"), self.option("evalue"), self.option("num_threads"))
        else:
            cmd = "{} -query {} -out {} -db {} -outfmt {} -evalue {} -num_threads {}"\
                .format(self.blastn_path, self.option("query_fa").prop["path"],
                        outfile_name, self.option("dbname_nsq"), self.option("outfmt"),
                        self.option("evalue"), self.option("num_threads"))
        if self.option("num_alignments"):
            cmd += " -num_alignments {}".format(self.option("num_alignments"))
        self.logger.info(cmd)
        command = self.add_command("blastn", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("blastn运行完成")
        else:
            self.set_error("blastn运行失败", code="34500601")

    def run(self):
        super(BlastNTool, self).run()
        self.run_blastn()
        self.end()
