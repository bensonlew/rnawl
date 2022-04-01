# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180521

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GenomeGrenameAgent(Agent):
    """
    带注释基因组的拆分
    """
    def __init__(self, parent):
        super(GenomeGrenameAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "string"},     # 下载基因组的fa文件
            {"name": "gff", "type": "string"},       # 下载基因组的gff文件
            {"name": "chromosomelist", "type": "string"}       # 更名文件
        ]
        self.add_option(options)
        self.step.add_steps('genomegrename')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.genomegrename.start()
        self.step.update()

    def step_end(self):
        self.step.genomegrename.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta"):
            raise OptionError("请设置fasta", code="34503301")
        if not self.option("gff"):
            raise OptionError("请设置 gff", code="34503302")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(GenomeGrenameAgent, self).end()


class GenomeGrenameTool(Tool):
    def __init__(self, config):
        super(GenomeGrenameTool, self).__init__(config)
        self.grename_path = self.config.PACKAGE_DIR + "/wgs/GRename.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def run_genomegrename(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{} {} -i {} -g {} -o {}"\
            .format(self.perl_path, self.grename_path, self.option("fasta"), self.option("gff"),
                    self.output_dir+"/ref")
        if self.option("chromosomelist"):
            cmd += " -f {} ".format(self.option("chromosomelist"))
        self.logger.info(cmd)
        self.logger.info("开始进行GenomeGrename")
        command = self.add_command("genomegrename", cmd).run()  # 必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("GenomeGrename完成！")
        else:
            self.set_error("GenomeGrename出错！", code="34503301")

    def run(self):
        super(GenomeGrenameTool, self).run()
        self.run_genomegrename()
        self.end()
