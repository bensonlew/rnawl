# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20181227

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class VcfConvertAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(VcfConvertAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_evolution.vcf"}
        ]
        self.add_option(options)
        self.step.add_steps('VcfConvert')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.VcfConvert.start()
        self.step.update()

    def step_end(self):
        self.step.VcfConvert.finish()
        self.step.update()

    def check_options(self):
        if not self.option('vcf_path'):
            raise OptionError('必须输入:vcf_path', code="35500403")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        super(VcfConvertAgent, self).end()


class VcfConvertTool(Tool):
    def __init__(self, config):
        super(VcfConvertTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl '
        self.vcf_convert = self.config.PACKAGE_DIR + '/noref_wgs/vcf-convert.pl'

    def run_vcf_convert(self):
        """
        """
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.vcf_convert, self.option("vcf_path").prop['path'],
                                         self.output_dir + "/pop.recode.vcf")
        self.logger.info("开始进行vcf_convert")
        command = self.add_command("vcf_convert", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("vcf_convert完成！")
        else:
            self.set_error("vcf_convert出错！", code="35500403")

    def run(self):
        super(VcfConvertTool, self).run()
        self.run_vcf_convert()
        self.end()
