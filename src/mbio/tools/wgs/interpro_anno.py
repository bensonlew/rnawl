# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20190225

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class InterproAnnoAgent(Agent):
    """
    参考基因组的PFAM数据库注释
    """
    def __init__(self, parent):
        super(InterproAnnoAgent, self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},
            {"name": "interpro_path", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('analysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.analysis.start()
        self.step.update()

    def step_end(self):
        self.step.analysis.finish()
        self.step.update()

    def check_options(self):
        if not self.option("fasta.list"):
            raise OptionError("请设置fasta.list")
        if not self.option("interpro_path"):
            raise OptionError("请设置interpro_path参数")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(InterproAnnoAgent, self).end()


class InterproAnnoTool(Tool):
    def __init__(self, config):
        super(InterproAnnoTool, self).__init__(config)
        # self.interpro_path = self.config.PACKAGE_DIR + "/wgs/interproanno.pl"
        self.interpro_path = self.config.PACKAGE_DIR + "/wgs/InterProanno_v2.pl"
        self.perl_path = 'miniconda2/bin/perl '

    def script_run(self):
        """
        perl interproanno.pl -falist fasta.list -anno 10.interPro -o /interPro.anno
        :return:
        """
        cmd = "{}{} -falist {} -anno {} -o {}"\
            .format(self.perl_path, self.interpro_path, self.option("fasta.list"), self.option("interpro_path"),
                    self.output_dir + "/test")
        self.logger.info(cmd)
        self.logger.info("开始进行interproAnno")
        command = self.add_command("interproanno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("interproanno完成！")
        else:
            self.set_error("interproanno出错！")

    def run(self):
        super(InterproAnnoTool, self).run()
        self.script_run()
        self.end()
