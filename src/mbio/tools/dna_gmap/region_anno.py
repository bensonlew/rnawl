
# -*- coding: utf-8 -*-
# __author__ = 'Zhao Binbin'
# modified 2018.06.22

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class RegionAnnoAgent(Agent):
    """
    关联区域注释
    """
    def __init__(self, parent):
        super(RegionAnnoAgent, self).__init__(parent)
        options = [
            {"name": "qtl_csv", "type": "infile", "format": "dna_gmap.lg"},
            {"name": "pop_final_vcf", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "pop_summary", "type": "infile", "format": "dna_gmap.lg"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("qtl_csv").is_set:
            raise OptionError("请设置性状的qtl.csv文件", code="34801601")
        if not self.option("pop_final_vcf").is_set:
            raise OptionError("请设置pop.final.vcf.gz文件", code="34801602")
        if not self.option("pop_summary").is_set:
            raise OptionError("请设置pop.summary文件", code="34801603")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(RegionAnnoAgent, self).end()


class RegionAnnoTool(Tool):
    def __init__(self, config):
        super(RegionAnnoTool, self).__init__(config)
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.region_qtl_extend_path = self.config.PACKAGE_DIR + "/dna_gmap/region-qtl.extend.pl"
        self.region_vcf_path = self.config.PACKAGE_DIR + "/dna_gmap/region-vcf.pl"
        self.region_gene_path = self.config.PACKAGE_DIR + "/dna_gmap/region-gene-new.anno.pl"

    def run_region_qtl_extend(self):
        self.trait = os.path.basename(self.option("qtl_csv").prop["path"]).split(".qtl.csv")[0]
        self.logger.info(self.trait)
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.region_qtl_extend_path, self.option("qtl_csv").prop["path"],
                                         self.output_dir + "/" + self.trait + ".qtl.result")
        command = self.add_command("qtl_extend", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本region_qtl_extend运行完成")
        else:
            self.set_error("脚本region_qtl_extend运行失败", code="34801601")
            self.set_error("脚本region_qtl_extend运行失败", code="34801606")

    def run_region_anno(self):
        cmd = "{} {} -i {} -o {} -r {}".format(self.perl_path, self.region_vcf_path, self.option("pop_final_vcf").prop["path"],
                                               self.output_dir + "/" + self.trait + ".vcf", self.output_dir + "/" + self.trait + ".qtl.result")
        command = self.add_command("region_anno", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本region_vcf运行完成")
        else:
            self.set_error("脚本region_vcf运行失败", code="34801602")

    def run_region_gene(self):
        cmd = "{} {} -i {} -a {} -o {}".format(self.perl_path, self.region_gene_path, self.output_dir + "/" + self.trait + ".qtl.result",
                                               self.option("pop_summary").prop["path"], self.output_dir + "/" + self.trait + ".gene")
        command = self.add_command("region_anno_db", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本region-gene-new.anno.pl运行完成")
        else:
            self.set_error("脚本region-gene-new.anno.pl运行失败", code="34801603")

    def run(self):
        super(RegionAnnoTool, self).run()
        self.run_region_qtl_extend()
        self.run_region_anno()
        self.run_region_gene()
        self.end()
