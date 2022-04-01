# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __modify__ = '2020/04/13'

from biocluster.workflow import Workflow
import os
from mbio.packages.bacgenome.common import fasta_cutoff


class GapFillWorkflow(Workflow):
    """
    
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GapFillWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "chr_fa", "type": "infile", "format": "sequence.fasta"},
            {'name': "cir_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "pre_log", "type": "infile", "format": "sequence.profile_table"},
            {"name": "genome_json", "type": "string"},
            {"name": "overlap", "type": "int"},
            {"name": "identity", "type": "float"},
            {"name": "select_seq", "type": "string"},
            {"name": "samp", "type": "string"},
            {"name": "gap", "type": "string"},
            {"name": "chr_len", "type": "string", "default": "1000"},
            {"name": "cir_len", "type": "string", "default": "1000"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.gap_fill = self.add_module('bacgenome.chr_gap_fill')
        self.seq_circle = self.add_tool('bacgenome.seq_circle')

    def run(self):
        if self.option('select_seq') in ['gap_fill']:
            self.run_gap_fill()
        elif self.option('select_seq') in ['assemble']:
            if self.option('gap'):
                self.run_seq_circle()
            else:
                self.run_gap_fill()
        super(GapFillWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.blast')
        self.logger.info("开始进行导表")
        seq_path = os.path.join(self.sheet.output, "result.fa")
        soft_ware_log = "\n" + "\n".join(self.gap_fill.gap_data_list)
        self.logger.info(soft_ware_log)
        api_path.add_gap_fill_detail(self.option("main_id"), self.option("samp"), self.gap_fill.option("log").prop["path"], seq_path, sw=soft_ware_log)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.gap_fill.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
            # ["assemble.stat.xls", "txt", '统计结果文件'],
            # ["assemble.summary.xls", "txt", "统计结果文件"],
            # ["all.fa", "", "基因组序列文件"]
        ])
        super(GapFillWorkflow, self).end()

    def run_seq_circle(self):
        """
        判断组装序列是环形还是线性
        :return:
        """
        fasta_cutoff(self.option('chr_fa').prop['path'], self.option("chr_len"), self.work_dir + "/" + "seq1.fa")
        opts = {
            "query": self.work_dir + "/" + "seq1.fa",
            "sample_name": self.option("samp")
        }
        self.seq_circle.set_options(opts)
        self.seq_circle.on("end", self.run_gap_fill)
        self.seq_circle.run()

    def run_gap_fill(self):
        fa = ''
        pre_log = ''
        if self.option('select_seq') in ['gap_fill']:
            fasta_cutoff(self.option('chr_fa').prop['path'], self.option("chr_len"), self.work_dir + "/" + "seq1.fa")
            fa = self.work_dir + "/" + "seq1.fa"
            pre_log = self.option("pre_log")
        elif self.option('select_seq') in ['assemble']:
            if self.option('gap') not in ['unicycler']:
                fa = self.seq_circle.option('fa_out')
                pre_log = self.seq_circle.option('cir_out')
            else:
                fasta_cutoff(self.option('chr_fa').prop['path'], self.option("chr_len"),
                             self.work_dir + "/" + "seq1.fa")
                fa = self.work_dir + "/" + "seq1.fa"
                pre_log = self.option("pre_log")
        fasta_cutoff(self.option('cir_fa').prop['path'], self.option("cir_len"), self.work_dir + "/" + "seq2.fa")
        opts = {
            "chr_fa": fa,
            "cir_fa": self.work_dir + "/" + "seq2.fa",
            "pre_log": pre_log,
            "genome_json": self.option("genome_json"),
            "overlap": self.option("overlap"),
            "identity": self.option("identity")
        }
        self.gap_fill.set_options(opts)
        self.gap_fill.on("end", self.set_db)
        self.gap_fill.run()
