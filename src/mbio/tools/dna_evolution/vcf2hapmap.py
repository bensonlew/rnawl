# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180822
# last modify 20180827

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class Vcf2hapmapAgent(Agent):
    """
    群体进化，vcf过滤后的vcf文件，用来生成map矩阵
    """
    def __init__(self, parent):
        super(Vcf2hapmapAgent, self).__init__(parent)
        options = [
            {"name": "recode_vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "ref", "type": "string", "default": "ref"},  # 20180827代码新加 cui
        ]
        self.add_option(options)
        self.step.add_steps('Vcf2hapmap')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Vcf2hapmap.start()
        self.step.update()

    def step_end(self):
        self.step.Vcf2hapmap.finish()
        self.step.update()

    def check_options(self):
        if not self.option("recode_vcf_path").is_set:
            raise OptionError("请设置recode_vcf_path")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(Vcf2hapmapAgent, self).end()


class Vcf2hapmapTool(Tool):
    def __init__(self, config):
        super(Vcf2hapmapTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.hapmap_path = self.config.PACKAGE_DIR + "/dna_evolution/vcf2hapmap.pl"

    def Vcf2hapmap(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{} {} -i {} -o {}"\
            .format(self.perl_path, self.hapmap_path, self.option("recode_vcf_path").prop['path'],
                    self.output_dir + "/pop.hapmap")
        cmd += " -ref {}".format(str(self.option('ref')))
        self.logger.info(cmd)
        self.logger.info("开始进行Vcf2hapmap")
        command = self.add_command("vcf2hapmap", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            # 多一个文件核查，否则检查pl脚本
            self.logger.info("Vcf2hapmap完成！")
        else:
            self.set_error("Vcf2hapmap出错！")
            raise Exception("Vcf2hapmap出错！")
        if not os.path.exists(self.output_dir + "/pop.hapmap"):
            self.set_error("pop.hapmap生成出错！")
        else:
            self.logger.info("pop.hapmap生成！")

    def run(self):
        super(Vcf2hapmapTool, self).run()
        self.Vcf2hapmap()
        self.end()
