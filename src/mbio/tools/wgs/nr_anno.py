# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180417

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class NrAnnoAgent(Agent):
    """
    参考基因组的NR数据库注释
    """
    def __init__(self,parent):
        super(NrAnnoAgent,self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},  # one col
            #{"name": "nr_path", "type": "infile","format": "wgs.nr_path"},  # name:.py的名字;type:infile/outfile;format:wgs/nr_path.py路径下的名字 || 如果有，就不用重新写check。
            {"name": "nr_path", "type": "string"},  # name:.py的名字;type:infile/outfile;format:wgs/nr_path.py路径下的名字 || 如果有，就不用重新写check。
        ]
        self.add_option(options)
        self.step.add_steps('NRanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.NRanno.start()
        self.step.update()

    def step_end(self):
        self.step.NRanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta.list"):
            raise OptionError("请设置fasta.list", code="34504101") # 必须有
        if not self.option("nr_path"):
            raise OptionError("请设置nr_path参数", code="34504102")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(NrAnnoAgent,self).end()
########

class NrAnnoTool(Tool):
    def __init__(self, config):
        super(NrAnnoTool, self).__init__(config)
        self.nr_path = self.config.PACKAGE_DIR + "/wgs/NRanno.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def NrAnno(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -d {} -o {}"\
            .format(self.perl_path, self.nr_path, self.option("fasta.list"), self.option("nr_path"),
                    self.output_dir + "/test")
        self.logger.info(cmd)
        self.logger.info("开始进行NrAnno")
        command = self.add_command("nranno", cmd).run()  # nranno必须小写，
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("NrAnno完成！")
        else:
            self.set_error("NrAnno出错！", code="34504101")
            self.set_error("NrAnno出错！", code="34504104")

    def run(self):
        super(NrAnnoTool, self).run()
        self.NrAnno()
        self.end()
