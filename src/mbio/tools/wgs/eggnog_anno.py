# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180423

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class EggnogAnnoAgent(Agent):
    """
    参考基因组的EGGNOG数据库注释
    """
    def __init__(self,parent):
        super(EggnogAnnoAgent,self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},
            {"name": "eggnog_path", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('EGGNOGanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.EGGNOGanno.start()
        self.step.update()

    def step_end(self):
        self.step.EGGNOGanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta.list"):
            raise OptionError("请设置fasta.list", code="34502301") # 必须有
        if not self.option("eggnog_path"):
            raise OptionError("请设置eggnog_path参数", code="34502302")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(EggnogAnnoAgent,self).end()
########

class EggnogAnnoTool(Tool):
    def __init__(self, config):
        super(EggnogAnnoTool, self).__init__(config)
        self.eggnog_path = self.config.PACKAGE_DIR + "/wgs/EGGanno.pl"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '

    def EggnogAnno(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -d {} -o {}"\
            .format(self.perl_path, self.eggnog_path, self.option("fasta.list"), self.option("eggnog_path"),
                    self.output_dir+"/test")
        self.logger.info(cmd)
        self.logger.info("开始进行EggnogAnno")
        command = self.add_command("eggnoganno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("EggnogAnno完成！")
        else:
            self.set_error("EggnogAnno出错！", code="34502301")
            self.set_error("EggnogAnno出错！", code="34502304")

    def run(self):
        super(EggnogAnnoTool, self).run()
        self.EggnogAnno()
        self.end()
