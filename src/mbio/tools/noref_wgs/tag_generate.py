# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20181227

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class TagGenerateAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(TagGenerateAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "dna_evolution.vcf"},
            {"name": "tsv_path", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('TagGenerate')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.TagGenerate.start()
        self.step.update()

    def step_end(self):
        self.step.TagGenerate.finish()
        self.step.update()

    def check_options(self):
        if not self.option('vcf_path'):
            raise OptionError('必须输入:vcf_path', code="35501305")
        if not self.option('tsv_path'):
            raise OptionError('必须输入:tsv_path', code="35501306")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(TagGenerateAgent, self).end()


class TagGenerateTool(Tool):
    def __init__(self, config):
        super(TagGenerateTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        self.tag_generate = self.config.PACKAGE_DIR + '/noref_wgs/tag-generate.pl'

    def run_tag_generate(self):
        """
        """
        cmd = "{} {} -vcf {} -catalog {} -out {}".format(self.perl_path, self.tag_generate,
                                                         self.option("vcf_path").prop['path'], self.option("tsv_path"),
                                                         self.output_dir + "/populations.tag")
        self.logger.info("开始进行tag_generate")
        command = self.add_command("tag_generate", cmd).run()  # 必须小写，
        self.wait()
        if command.return_code == 0:
            self.logger.info("tag_generate完成！")
        else:
            self.set_error("tag_generate出错！", code="35501303")

    def run(self):
        super(TagGenerateTool, self).run()
        self.run_tag_generate()
        self.end()
