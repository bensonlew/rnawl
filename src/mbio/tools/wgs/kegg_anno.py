# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180423

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class KeggAnnoAgent(Agent):
    """
    参考基因组的KEGG数据库注释
    """
    def __init__(self,parent):
        super(KeggAnnoAgent,self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},
            {"name": "kegg_path", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('KEGGanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.KEGGanno.start()
        self.step.update()

    def step_end(self):
        self.step.KEGGanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta.list"):
            raise OptionError("请设置fasta.list", code="34503701") # 必须有
        if not self.option("kegg_path"):
            raise OptionError("请设置kegg_path参数", code="34503702")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(KeggAnnoAgent,self).end()
########

class KeggAnnoTool(Tool):
    def __init__(self, config):
        super(KeggAnnoTool, self).__init__(config)
        self.kegg_path = self.config.PACKAGE_DIR + "/wgs/KEGGanno.pl"
        self.perl_path = 'miniconda2/bin/perl '

    def KeggAnno(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -d {} -o {}"\
            .format(self.perl_path, self.kegg_path, self.option("fasta.list"), self.option("kegg_path"),
                    self.output_dir + "/test")
        self.logger.info(cmd)
        self.logger.info("开始进行KeggAnno")
        command = self.add_command("kegganno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("KeggAnno完成！")
        else:
            self.set_error("KeggAnno出错！", code="34503701")
            self.set_error("KeggAnno出错！", code="34503704")

    def run(self):
        super(KeggAnnoTool, self).run()
        self.KeggAnno()
        self.end()
