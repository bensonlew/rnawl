# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20181226

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import shutil


class GenotypeAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GenotypeAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "string"},   # 输入输出路径
            {"name": "gro_list", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('Genotype')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Genotype.start()
        self.step.update()

    def step_end(self):
        self.step.Genotype.finish()
        self.step.update()

    def check_options(self):
        if not self.option('input_file'):
            raise OptionError('必须输入:input_file', code="300105")
        if not self.option('gro_list'):
            raise OptionError('必须输入:gro_list', code="300106")

    def set_resource(self):
        """
        运行所需资源
        vcftools一般不耗费内存资源，这里暂时不设置动态内存
        """
        self._cpu = 20
        self._memory = "180G"

    def end(self):
        super(GenotypeAgent, self).end()


class GenotypeTool(Tool):
    def __init__(self, config):
        super(GenotypeTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.tsv2bam_path = "bioinfo/noRefWGS/tsv2bam"
        self.gstacks_path = "bioinfo/noRefWGS/gstacks"
        self.populations_path = "bioinfo/noRefWGS/populations"

    def run_tsv2bam(self):
        """
        """
        cmd = "{} -P {} -M {} -t 16".format(self.tsv2bam_path, self.option("input_file"), self.option("gro_list"))
        self.logger.info("开始进行Tsv2Bam")
        command = self.add_command("tsv2bam", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("Tsv2Bam完成！")
        else:
            self.set_error("Tsv2Bam出错！", code="300107")

    def run_gstacks(self):
        """
        """
        cmd = "{} -P {} -M {} -t 16".format(self.gstacks_path, self.option("input_file"), self.option("gro_list"))
        self.logger.info("开始进行gstacks")
        command = self.add_command("gstacks", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("gstacks完成！")
        else:
            self.set_error("gstacks出错！", code="300108")

    def run_populations(self):
        """
        {} -P {} -M {} -t 16 -O {} --vcf --fasta_loci --hwe --fstats -k
        modified by hd 为了保持与线下代码一致，删除了--hwe --fstats -k
        modified by hd@20200525 再次与黄总监线下流程一致，增加--hwe --fstats -k，同时stacks升级后，
        --fasta_loci改成了--fasta-loci
        """
        cmd = "{} -P {} -M {} -t 16 -O {} --vcf --fasta-loci --hwe --fstats".format(
            self.populations_path, self.option("input_file"), self.option("gro_list"), self.option("input_file"))
        self.logger.info("开始进行populations")
        command = self.add_command("populations", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("populations完成！")
        else:
            self.set_error("populations出错！", code="300109")

    def cp_outfile(self):
        vcf_path = os.path.join(self.option("input_file"), "populations.snps.vcf")
        tsv_path = os.path.join(self.option("input_file"), "catalog.tags.tsv.gz")
        output_vcf_path = os.path.join(self.output_dir, "populations.snps.vcf")
        output_tsv_path = os.path.join(self.output_dir, "catalog.tags.tsv.gz")
        if os.path.exists(output_vcf_path):
            os.remove(output_vcf_path)
        if os.path.exists(output_tsv_path):
            os.remove(output_tsv_path)
        os.link(vcf_path, output_vcf_path)
        os.link(tsv_path, output_tsv_path)

    def run(self):
        super(GenotypeTool, self).run()
        self.run_tsv2bam()
        self.run_gstacks()
        self.run_populations()
        self.cp_outfile()
        self.end()
