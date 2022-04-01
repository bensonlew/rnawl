# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180611

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GetgenefastaAgent(Agent):
    """
    提取ref.fa里面的gene.fa文件
    """
    def __init__(self,parent):
        super(GetgenefastaAgent,self).__init__(parent)
        options = [
            {"name": "reffa", "type": "string"}, # 更名后的reffa
            {"name": "refgff", "type": "string"}, # 更名后的refgff
        ]
        self.add_option(options)
        self.step.add_steps('Getgenefasta')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Getgenefasta.start()
        self.step.update()

    def step_end(self):
        self.step.Getgenefasta.finish()
        self.step.update()

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置reffa参数", code="24500801")
        if not self.option("refgff"):
            raise OptionError("请设置refgff参数", code="24500802")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(GetgenefastaAgent,self).end()
########

class GetgenefastaTool(Tool):
    def __init__(self, config):
        super(GetgenefastaTool, self).__init__(config)
        self.gene_path = self.config.PACKAGE_DIR + "/wgs/getGeneFasta.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def Getgenefasta(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -g {} -o {}"\
            .format(self.perl_path, self.gene_path, self.option("reffa"), self.option("refgff"),
                    self.output_dir + "/ref.gene.fa")
        self.logger.info(cmd)
        self.logger.info("开始进行Getgenefasta")
        command = self.add_command("getgenefasta", cmd).run() 
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Getgenefasta完成！")
        else:
            self.set_error("Getgenefasta出错！", code="24500801")
            self.set_error("Getgenefasta出错！", code="24500804")

    def run(self):
        super(GetgenefastaTool, self).run()
        self.Getgenefasta()
        self.end()
