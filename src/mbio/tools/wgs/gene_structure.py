# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.07

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GeneStructureAgent(Agent):
    """
    基因详情页-基因结构
    """
    def __init__(self, parent):
        super(GeneStructureAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "string"},
            {"name": "ref_gff", "type": "string"},
            {"name": "pop_final_vcf", "type": "infile", 'format': 'bsa.vcf'},
            {"name": "gene_id", "type": "string"},
            {"name": "chrom", "type": "string"},  # 基因id所在的染色体位置
            {"name": "location", "type": "string"},  # 基因位置
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ref_fa"):
            raise OptionError("请设置ref_fa文件", code="34503201")
        if not self.option("ref_gff"):
            raise OptionError("请设置ref_gff文件", code="34503202")
        if not self.option("pop_final_vcf"):
            raise OptionError("请设置pop_final_vcf文件", code="34503203")
        if not self.option("gene_id"):
            raise OptionError("请设置gene_id", code="34503204")
        if not self.option("location"):
            raise OptionError("请设置基因位置location", code="34503205")

    def set_resource(self):
        self._cpu = 3
        size = os.path.getsize(self.option("pop_final_vcf").prop['path'])
        size = size / 1024 / 1024 / 1024
        if size < 10:
            self._memory = "15G"
        elif size < 20:
            self._memory = "25G"
        elif size < 30:
            self._memory = "35G"
        else:
            self._memory = "50G"

    def end(self):
        super(GeneStructureAgent, self).end()


class GeneStructureTool(Tool):
    def __init__(self, config):
        super(GeneStructureTool, self).__init__(config)
        self.python = "miniconda2/bin/python"
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.exon_seq = self.config.PACKAGE_DIR + "/wgs/gene_exon_seq.pl"
        self.gene_alt = self.config.PACKAGE_DIR + "/wgs/gene_alt.py"
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome"
        self.ref_fa = os.path.join(self.json_path, self.option("ref_fa"))
        self.ref_gff = os.path.join(self.json_path, self.option("ref_gff"))

    def run_exon_seq(self):
        """
        运行gene_exon_seq.pl，得到ref的序列信息
        """
        cmd = "{} {} -ref {} -gff {}".format(self.perl_path, self.exon_seq, self.ref_fa, self.ref_gff)
        cmd += " -gene {} -out {}".format(self.option("gene_id"), self.output_dir)
        command = self.add_command("gene_exon_seq", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gene_exon_seq.pl运行成功")
        else:
            self.set_error("gene_exon_seq.pl运行失败", code="34503201")

    def run_gene_alt(self):
        """
        运行gene_alt.py，得到alt的信息
        """
        start = self.option("location").split("-")[0]
        end = self.option("location").split("-")[1]
        cmd = "{} {} -gene {} -chr {}".format(self.python, self.gene_alt, self.option("gene_id"), self.option("chrom"))
        cmd += " -start {} -end {} -vcf {}".format(start, end, self.option("pop_final_vcf").prop['path'])
        cmd += " -out {}".format(os.path.join(self.output_dir, self.option("gene_id") + ".xls"))
        command = self.add_command("gene_alt", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gene_alt.py运行成功")
        else:
            self.set_error("gene_alt.py运行失败", code="34503202")

    def set_db(self):
        """
        将结果保存到mongo
        """
        seq_api = self.api.api("wgs.gene_detail")
        if self.option("project_type"):
            seq_api._project_type = self.option("project_type")
        structure_id = self.option("main_id")
        exon_path = os.path.join(self.output_dir, self.option("gene_id") + ".exon.fa")
        mirna_path = os.path.join(self.output_dir, self.option("gene_id") + ".mirna.xls")
        alt_path = os.path.join(self.output_dir, self.option("gene_id") + ".xls")
        gene_id = self.option("gene_id")
        start = int(self.option("location").split("-")[0])
        end = int(self.option("location").split("-")[1])
        seq_api.add_sg_gene_structure_detail(structure_id, gene_id, start, end, exon_path, alt_path, mirna_path)

    def run(self):
        super(GeneStructureTool, self).run()
        self.run_exon_seq()
        self.run_gene_alt()
        self.set_db()
        self.end()
