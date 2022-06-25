# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20190225

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class PfamAnnoAgent(Agent):
    """
    参考基因组的PFAM数据库注释
    """
    def __init__(self, parent):
        super(PfamAnnoAgent, self).__init__(parent)
        options = [
            {"name": "fasta.list", "type": "string"},
            {"name": "pfam_path", "type": "string"},
            {"name": "anno_method", "type": "string", "default": "interpro_scan"}  # interpro_scan or pfam_scan
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
        if not self.option("pfam_path"):
            raise OptionError("请设置pfam_path参数")
        if self.option("anno_method") not in ['pfam_scan', 'interpro_scan']:
            raise OptionError("anno method is not right, must be pfam_scan or interpro_scan")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(PfamAnnoAgent, self).end()


class PfamAnnoTool(Tool):
    def __init__(self, config):
        super(PfamAnnoTool, self).__init__(config)
        if self.option("anno_method") == 'pfam_scan':
            self.pfam_path = self.config.PACKAGE_DIR + "/wgs/pfamanno.pl"
        else:
            self.pfam_path = self.config.PACKAGE_DIR + "/wgs/Pfamanno_v2.pl"
        self.perl_path = 'miniconda2/bin/perl '

    def script_run(self):
        """
        perl pfamanno.pl -falist fasta.list -anno 09.Pfam -o pfam.anno
        :return:
        """
        cmd = "{}{} -falist {} -anno {} -o {}"\
            .format(self.perl_path, self.pfam_path, self.option("fasta.list"), self.option("pfam_path"),
                    self.output_dir + "/test")
        self.logger.info(cmd)
        self.logger.info("开始进行pfamAnno")
        command = self.add_command("pfamanno", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pfamAnno完成！")
        else:
            self.set_error("pfamAnno出错！")

    def run(self):
        super(PfamAnnoTool, self).run()
        self.script_run()
        self.end()
