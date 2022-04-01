# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180822

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class Vcf2tritAgent(Agent):
    """
    GWAS流程的trit.pl,输入trait文件，输出一个文件夹的性状文件。
    """
    def __init__(self, parent):
        super(Vcf2tritAgent, self).__init__(parent)
        options = [
            {"name": "upload_trait_path", "type": "infile", "format": "dna_gmap.trait"},
        ]
        self.add_option(options)
        self.step.add_steps('Vcf2trit')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Vcf2trit.start()
        self.step.update()

    def step_end(self):
        self.step.Vcf2trit.finish()
        self.step.update()

    def check_options(self):
        if not self.option("upload_trait_path").is_set:
            raise OptionError("请上传upload_trait_path")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(Vcf2tritAgent, self).end()


class Vcf2tritTool(Tool):
    def __init__(self, config):
        super(Vcf2tritTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.trait_path = self.config.PACKAGE_DIR + "/dna_evolution/trit.pl"

    def Vcf2trit(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{} {} -trt {} -out {}"\
            .format(self.perl_path, self.trait_path, self.option("upload_trait_path").prop['path'],
                    self.output_dir + "/")
        self.logger.info(cmd)
        self.logger.info("开始进行Vcf2trit")
        command = self.add_command("vcf2trit", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Vcf2trit完成！")
        else:
            self.set_error("Vcf2trit出错！")
            raise Exception("Vcf2trit出错！")

    def run(self):
        super(Vcf2tritTool, self).run()
        self.Vcf2trit()
        self.end()
