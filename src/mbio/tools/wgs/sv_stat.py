# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.10

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SvStatAgent(Agent):
    """
    sv 变异类型统计,变异长度统计
    """
    def __init__(self, parent):
        super(SvStatAgent, self).__init__(parent)
        options = [
            {"name": "sv_anno", "type": "string"},  # sv注释文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sv_anno"):
            raise OptionError("请设置sv注释文件", code="34507001")
        if not os.path.exists(self.option("sv_anno")):
            raise OptionError("sv注释文件:%s不存在，请检查", variables=(self.option("sv_anno")), code="34507002")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SvStatAgent, self).end()


class SvStatTool(Tool):
    def __init__(self, config):
        super(SvStatTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.sv_stat = self.config.PACKAGE_DIR + "/wgs/sv_stat.pl"

    def run_sv_stat(self):
        """
        sv_stat.pl sv变异类型和长度统计
        """
        self.sample_name = os.path.basename(self.option("sv_anno")).split(".sv.anno")[0]
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.sv_stat, self.option("sv_anno"), self.work_dir + "/" + self.sample_name)
        command = self.add_command("sv_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sv_stat.pl运行完成")
        else:
            self.set_error("sv_stat.pl运行失败", code="34507001")

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.work_dir + "/" + self.sample_name + ".sv.stat", self.output_dir + "/" + self.sample_name + ".sv.stat.xls")
        os.link(self.work_dir + "/" + self.sample_name + ".sv.length", self.output_dir + "/" + self.sample_name + ".sv.length.xls")

    def run(self):
        super(SvStatTool, self).run()
        self.run_sv_stat()
        self.set_output()
        self.end()
