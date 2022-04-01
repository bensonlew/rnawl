# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2020/03/31'

from biocluster.workflow import Workflow
from mbio.packages.bacgenome.common import fasta_cutoff,get_scaffoldlist
import os


class ScaffoldBaseCoverageWorkflow(Workflow):
    """
    用于扫描图的scaffolds base coverage分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ScaffoldBaseCoverageWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},
            {"name": "samp", "type": "string"},
            {"name": "min_len", "type": "string", "default":"1000"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.bowtie2 = self.add_tool("align.bowtie2")
        self.convert = self.add_tool('metagbin.convert_format')
        self.coverage = self.add_tool("bacgenome.scaf_mean_cov")
        self.file_path = self._sheet.output

    def check_options(self):
        """
        检查参数
        :return:
        """
        self.prefix = os.path.basename(self.option("read1").prop["path"]).split(".")[0]
        return True

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_bowtie()
        super(ScaffoldBaseCoverageWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('bac_assem.assemble')
        self.logger.info("开始进行导表")
        sample = self.option("samp")
        depth_detail = os.path.join(self.coverage.output_dir, "cov_mean.txt")
        listseqid = get_scaffoldlist(self.work_dir + "/" + "seq.fa")
        api_path.add_draft_seq_cov(self.option("main_id"), depth_detail, sample, listseqid)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.coverage.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果目录"],
        ])
        super(ScaffoldBaseCoverageWorkflow, self).end()

    def run_bowtie(self):
        fasta_cutoff(self.option('fa').prop['path'], self.option("min_len"), self.work_dir + "/" + "seq.fa")
        self.bowtie2.set_options({
            "fastq1": self.option("read1"),
            "fastq2": self.option("read2"),
            "ref_fasta": self.work_dir + "/" + "seq.fa"
        })
        self.bowtie2.on("end", self.run_convert)
        self.bowtie2.run()

    def run_convert(self):
        self.convert.set_options({
            "sam": self.bowtie2.output_dir + "/" + self.prefix + ".pair.sam",
            "analysis": "metagbin"
        })
        self.convert.on("end", self.run_coverage)
        self.convert.run()

    def run_coverage(self):
        self.coverage.set_options({
            "bam_file": self.convert.output_dir + "/" + self.prefix + "_sorted.bam"
        })
        self.coverage.on("end", self.set_db)
        self.coverage.run()
