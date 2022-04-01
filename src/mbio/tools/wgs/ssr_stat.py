# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SsrStatAgent(Agent):
    """
    软件：SSR.type.stat.pl
    """
    def __init__(self, parent):
        super(SsrStatAgent, self).__init__(parent)
        options = [
            {"name": "misa_file", "type": "string"},  # scafSeq.misa文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("misa_file"):
            raise OptionError("请设置misa_file文件", code="34506601")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SsrStatAgent, self).end()


class SsrStatTool(Tool):
    def __init__(self, config):
        super(SsrStatTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.ssr_type_stat = self.config.PACKAGE_DIR + "/wgs/SSR.type.stat.pl"

    def run_ssr_type_stat(self):
        """
        SSR.type.stat.pl
        """
        ssr_stat = self.output_dir + "/" + os.path.basename(self.option("misa_file")).split(".")[0] + ".ssr.stat"
        cmd = "{} {}  -m {} -o {}".format(self.perl_path, self.ssr_type_stat, self.option("misa_file"), ssr_stat)
        command = self.add_command("ssr_type_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr_type_stat完成")
        else:
            self.set_error("ssr_type_stat失败", code="34506601")

    def run(self):
        super(SsrStatTool, self).run()
        self.run_ssr_type_stat()
        self.end()
