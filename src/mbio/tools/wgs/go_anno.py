# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180423

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GoAnnoAgent(Agent):
    """
    参考基因组的GO数据库注释
    """
    def __init__(self,parent):
        super(GoAnnoAgent,self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},
            {"name": "go_path", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('GOanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.GOanno.start()
        self.step.update()

    def step_end(self):
        self.step.GOanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta.list"):
            raise OptionError("请设置fasta.list", code="34503401") # 必须有
        if not self.option("go_path"):
            raise OptionError("请设置go_path参数", code="34503402")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(GoAnnoAgent,self).end()
########

class GoAnnoTool(Tool):
    def __init__(self, config):
        super(GoAnnoTool, self).__init__(config)
        self.go_path = self.config.PACKAGE_DIR + "/wgs/GOanno.pl"
        self.perl_path = 'miniconda2/bin/perl '

    def GoAnno(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -d {} -o {}"\
            .format(self.perl_path, self.go_path, self.option("fasta.list"), self.option("go_path"),
                    self.output_dir+ "/test")
        self.logger.info(cmd)
        self.logger.info("开始进行GoAnno")
        command = self.add_command("goanno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("GoAnno完成！")
        else:
            self.set_error("GoAnno出错！", code="34503401")
            self.set_error("GoAnno出错！", code="34503404")

    def run(self):
        super(GoAnnoTool, self).run()
        self.GoAnno()
        self.end()
