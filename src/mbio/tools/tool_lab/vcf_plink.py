# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re

class VcfPlinkAgent(Agent):
    def __init__(self, parent):
        super(VcfPlinkAgent, self).__init__(parent)
        options = [
            {"name":"vcf_file","type":"infile","format":"dna_gmap.vcf"},
            {'name': 'chrom_map', "type": 'infile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vcf_file"):
            raise OptionError("请设置VCF文件")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 8
        size = os.path.getsize(self.option("vcf_file").prop['path'])
        size = size / 1024 / 1024 / 1024
        if size < 50:
            self._memory = "40G"
        elif size < 80:
            self._memory = "70G"
        elif size < 100:
            self._memory = "90G"
        else:
            self._memory = "120G"

    def end(self):
        super(VcfPlinkAgent, self).end()

class VcfPlinkTool(Tool):
    def __init__(self, config):
        super(VcfPlinkTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.vcftools_path = "bioinfo/dna_evolution/vcftools"
        self.plink_path = "bioinfo/dna_evolution/plink/plink"
    
    def run_vcf_plink(self):
        cmd = "{} --vcf {}".format(self.vcftools_path, self.option('vcf_file').prop['path'])
        if self.option('chrom_map'):
            cmd += " --plink --chrom-map {} --out {}".format(self.option('chrom_map').prop['path'],self.output_dir + "/pop")
        else:
            cmd += " --plink --out {}".format(self.output_dir + "/pop")
        self.logger.info(cmd)
        self.logger.info("开始进行VcftoolsPlink")
        command = self.add_command("vcftoolsplink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("VcftoolsPlink完成！")
        else:
            self.set_error("VcftoolsPlink出错！")

    def run_plink(self):
        cmd = "{}".format(self.plink_path)
        cmd += " --file {} ".format(self.output_dir + "/pop")
        cmd += " --make-bed"
        cmd += " --out {}".format(self.output_dir + "/pop")
        self.logger.info(cmd)
        self.logger.info("开始进行Plink")
        command = self.add_command("plink", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Plink完成！")
        else:
            self.set_error("Plink出错！")

    def run(self):
        super(VcfPlinkTool, self).run()
        self.run_vcf_plink()
        self.run_plink()
        self.end()