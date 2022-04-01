# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180427

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SplitFastawgsAgent(Agent):
    """
    带注释基因组的拆分
    """
    def __init__(self, parent):
        super(SplitFastawgsAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('Splitfastawgs')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Splitfastawgs.start()
        self.step.update()

    def step_end(self):
        self.step.Splitfastawgs.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta"):
            raise OptionError("请设置fasta", code="34506101")  # 必须有

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(SplitFastawgsAgent, self).end()


class SplitFastawgsTool(Tool):
    def __init__(self, config):
        super(SplitFastawgsTool, self).__init__(config)
        self.split_path = self.config.PACKAGE_DIR + "/wgs/splitFasta.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def splitfastawgs(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -n 20 -o {}"\
            .format(self.perl_path, self.split_path, self.option("fasta"),
                    self.output_dir+"/")
        self.logger.info(cmd)
        self.logger.info("开始进行SplitFastawgs")
        command = self.add_command("splitfastawgs", cmd).run()  # splitfastawgs必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("SplitFastawgs完成！")
        else:
            self.set_error("SplitFastawgs出错！", code="34506101")

    def run(self):
        super(SplitFastawgsTool, self).run()
        self.splitfastawgs()
        self.end()
