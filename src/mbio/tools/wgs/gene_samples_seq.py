# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.04

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GeneSamplesSeqAgent(Agent):
    """
    基因详情页-样本序列
    """
    def __init__(self, parent):
        super(GeneSamplesSeqAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "string"},
            {"name": "ref_gff", "type": "string"},
            {"name": "specimen_ids", "type": "string"},
            {"name": "pop_final_vcf", "type": "infile", 'format': 'bsa.vcf'},  # pop_final_vcf
            {"name": "gene_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa文件", code="34503101")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref_gff文件", code="34503102")
        if not self.option("specimen_ids"):
            raise OptionError("请设置specimen_ids", code="34503103")
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="34503104")
        if not self.option("gene_id"):
            raise OptionError("请设置gene_id", code="34503105")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(GeneSamplesSeqAgent, self).end()


class GeneSamplesSeqTool(Tool):
    def __init__(self, config):
        super(GeneSamplesSeqTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.gene_samples_seq = self.config.PACKAGE_DIR + "/wgs/gene_samples_seq.pl"
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome"  # 参考组配置文件
        self.ref_fa = os.path.join(self.json_path, self.option("ref_fa"))
        self.ref_gff = os.path.join(self.json_path, self.option("ref_gff"))

    def run_gene_samples_seq(self):
        """
        运行gene_samples_seq.pl脚本，得到样本列表对应的gene_id序列信息
        """
        with open(self.work_dir + "/sample.list", "w") as w:
            for s in self.option("specimen_ids").split(";"):
                w.write(s + "\n")
        cmd = "{} {} -ref {}".format(self.perl_path, self.gene_samples_seq, self.ref_fa)
        cmd += " -gff {} -vcf {}".format(self.ref_gff, self.option("pop_final_vcf").prop['path'])
        cmd += " -list {} -gene {} -out {}".format(self.work_dir + "/sample.list", self.option("gene_id"),
                                                   self.output_dir)
        command = self.add_command("get_samples_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gene_samples_seq.pl 运行完成")
        else:
            self.set_error("gene_samples_seq.pl 运行失败", code="34503101")

    def set_db(self):
        """
        保存结果到mongo
        """
        seq_api = self.api.api("wgs.gene_detail")
        if self.option("project_type"):
            seq_api._project_type = self.option("project_type")
        for f in os.listdir(self.output_dir):
            gene_path = os.path.join(self.output_dir, f)
        seq_api.add_sg_gene_sample_seq_detail(seq_id=self.option("main_id"), gene_id=self.option("gene_id"),
                                              gene_path=gene_path)

    def run(self):
        super(GeneSamplesSeqTool, self).run()
        self.run_gene_samples_seq()
        if self.option("main_id"):
            self.set_db()
        self.end()
