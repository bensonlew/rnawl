# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/5/9'

from biocluster.workflow import Workflow
import os
from mbio.packages.metagenomic.common import link_file,link_dir


class CompleteStatWorkflow(Workflow):
    """
    
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompleteStatWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "fq1", "type": "infile", "format": "sequence.fastq"},
            {"name": "fq2", "type": "infile", "format": "sequence.fastq"},
            {"name": "genome_id", "type": "string"},
            {"name": "samp", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.get_seq = self.add_tool("bacgenome.complete_get_seq")
        self.correct = self.add_module("bacgenome.pilon_correct")
        self.stat = self.add_tool('bacgenome.complete_asse_stat2')

    def run(self):
        self.run_get_seq()
        super(CompleteStatWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.blast')
        self.logger.info("开始进行导表")
        seq_path = os.path.join(self.sheet.output, "all.fa")
        table1 = os.path.join(self.stat.output_dir, "assemble.stat.xls")
        table2 = os.path.join(self.stat.output_dir, "assemble.summary.xls")
        api_path.add_complete_stat_detail(self.option("main_id"), self.option("genome_id"), self.option("samp"), seq_path,
                                  table1, table2)
        self.end()

    def end(self):
        link_dir(self.stat.output_dir, self.output_dir)
        link_file(self.correct.option("scaffold").prop["path"], self.output_dir + "/all.fa")
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
            ["assemble.stat.xls", "txt", '统计结果文件'],
            ["assemble.summary.xls", "txt", "统计结果文件"],
            ["all.fa", "", "基因组序列文件"],
            ["seq_dir", "", "基因组序列目录"],
        ])
        super(CompleteStatWorkflow, self).end()

    def run_get_seq(self):
        opts = {
            "fa": self.option("fa"),
            "genome_id": self.option("genome_id")
        }
        self.get_seq.set_options(opts)
        self.get_seq.on("end", self.run_correct)
        self.get_seq.run()

    def run_correct(self):
        opts = {
            "seq_scaf": self.get_seq.option("out_fa"),
            "fq1": self.option("fq1"),
            "fq2": self.option("fq2"),
            "is_major_result": True,
            "sample_name": self.option("samp")
        }
        self.correct.set_options(opts)
        self.correct.on("end", self.run_stat)
        self.correct.run()

    def run_stat(self):
        opts = {
            "fa": self.correct.option("scaffold"),
            "sample_name": self.option("samp")
        }
        self.stat.set_options(opts)
        self.stat.on("end", self.set_db)
        self.stat.run()
