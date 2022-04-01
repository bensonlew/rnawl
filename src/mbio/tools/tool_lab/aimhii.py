# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.03.12

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,shutil


class AimhiiAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(AimhiiAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # ref.fa
            {"name": "insert_fa", "type": "infile", "format": "sequence.fasta"},  # 外源片段序列
            {"name": "clean_fq_1", "type": "infile", "format": "sequence.fastq"},  # 279_1.clean.1.fastq.gz
            {"name": "clean_fq_2", "type": "infile", "format": "sequence.fastq"},  # 279_1.clean.2.fastq.gz
            {"name": "sample", "type": "string"},  # 样本名称
            {"name": "num", "type": "string", "default": 8},  # threads
            {"name": "result_csv", "type": "outfile", "format": "wgs_v2.bcf"}  # 279_1.result.csv
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa").is_set:
            raise OptionError("请设置ref_fa")
        if not self.option("insert_fa").is_set:
            raise OptionError("请设置insert_fa")
        if not self.option("clean_fq_1").is_set:
            raise OptionError("请设置clean_fq_1")
        if not self.option("sample"):
            raise OptionError("请设置sample")

    def set_resource(self):
        self._cpu = 5
        self._memory = "100G"

    def end(self):
        super(AimhiiAgent, self).end()


class AimhiiTool(Tool):
    def __init__(self, config):
        super(AimhiiTool, self).__init__(config)
        #self.aimhii = "/bioinfo/tool_lab/aimhii/aimhii_venv/bin/aimhii"
        self.aimhii = "program/Python/bin/aimhii"
        self.samtools = "program/Python/bin/samtools"
        self.illumina_adapter1 = self.config.SOFTWARE_DIR + "/database/wgs_v2/illumina_adapter1.fasta"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/wgs_v2')

    def run_aimhii(self):
        """

        """
        if os.path.exists(self.work_dir + "/tmp"):
            shutil.rmtree(self.work_dir + "/tmp")
        os.mkdir(self.work_dir + "/tmp")
        outfile = os.path.join(self.output_dir, (self.option("sample") + ".result.csv"))
        plot = os.path.join(self.output_dir, (self.option("sample") + ".readplot.pdf"))
        if self.option("clean_fq_2").is_set:
            cmd = "{} {} {} {} {} {} --outfile {} --plot {} --threads {} --tmpdir {} "\
                .format(self.aimhii, self.option("ref_fa").prop["path"], self.option("insert_fa").prop["path"],
                        self.illumina_adapter1, self.option("clean_fq_1").prop["path"],
                        self.option("clean_fq_2").prop["path"], outfile, plot, self.option("num"), self.work_dir + "/tmp")
        else:
            cmd = "{} {} {} {} {} --outfile {} --plot {} --threads {} --tmpdir {} "\
                .format(self.aimhii, self.option("ref_fa").prop["path"], self.option("insert_fa").prop["path"],
                         self.illumina_adapter1, self.option("clean_fq_1").prop["path"],
                        outfile, plot, self.option("num") ,self.work_dir + "/tmp")
        command = self.add_command("aimhii", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("aimhii运行成功")
        else:
            self.set_error("aimhii运行失败")
        self.option("result_csv", outfile)

    def run(self):
        super(AimhiiTool, self).run()
        self.run_aimhii()
        self.end()
