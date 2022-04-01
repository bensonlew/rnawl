# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180515

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SnpeffGffcheckAgent(Agent):
    """
    参考基因组snpeff check ref.gff生成genes.gff
    """
    def __init__(self,parent):
        super(SnpeffGffcheckAgent,self).__init__(parent)
        options = [
            {"name": "refgff", "type": "string"},  # one col
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("refgff"):
            raise OptionError("请输入ref.gff", code="34505801") # 必须有

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(SnpeffGffcheckAgent,self).end()
########

class SnpeffGffcheckTool(Tool):
    def __init__(self, config):
        super(SnpeffGffcheckTool, self).__init__(config)
        self.gffcheck_path = self.config.PACKAGE_DIR + "/wgs/gffcheck.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def SnpeffGffcheck(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{} {} -i {}  -o {}"\
            .format(self.perl_path, self.gffcheck_path, self.option("refgff"),
                    self.output_dir + "/genes.gff")
        self.logger.info(cmd)
        self.logger.info("开始进行SnpeffGffcheck")
        command = self.add_command("snpeffgffcheck", cmd).run()  # nranno必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("SnpeffGffcheck完成！")
        else:
            self.set_error("SnpeffGffcheck出错！", code="34505801")
            self.set_error("SnpeffGffcheck出错！", code="34505804")

    def run(self):
        super(SnpeffGffcheckTool, self).run()
        self.SnpeffGffcheck()
        self.end()
