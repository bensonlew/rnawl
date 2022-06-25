# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last modify 20180423

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class UniprotAnnoAgent(Agent):
    """
    参考基因组的UNIPROT数据库注释
    """
    def __init__(self,parent):
        super(UniprotAnnoAgent,self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},
            {"name": "uniprot_path", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('UNIPROTanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.UNIPROTanno.start()
        self.step.update()

    def step_end(self):
        self.step.UNIPROTanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta.list"):
            raise OptionError("请设置fasta.list", code="34507201") # 必须有
        if not self.option("uniprot_path"):
            raise OptionError("请设置uniprot_path参数", code="34507202")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(UniprotAnnoAgent,self).end()
########

class UniprotAnnoTool(Tool):
    def __init__(self, config):
        super(UniprotAnnoTool, self).__init__(config)
        self.uniprot_path = self.config.PACKAGE_DIR + "/wgs/UNIPROTanno.pl"
        self.perl_path = 'miniconda2/bin/perl '

    def UniprotAnno(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -d {} -o {}"\
            .format(self.perl_path, self.uniprot_path, self.option("fasta.list"), self.option("uniprot_path"),
                    self.output_dir+"/test")
        self.logger.info(cmd)
        self.logger.info("开始进行UniprotAnno")
        command = self.add_command("uniprotanno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("UniprotAnno完成！")
        else:
            self.set_error("UniprotAnno出错！", code="34507201")
            self.set_error("UniprotAnno出错！", code="34507204")

    def run(self):
        super(UniprotAnnoTool, self).run()
        self.UniprotAnno()
        self.end()
