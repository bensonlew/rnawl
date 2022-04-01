# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.12.11
# recode: modified by qinwen 2021.07.02

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime

starttime = datetime.datetime.now()


class LdDrawAgent(Agent):
    """
    连锁不平衡接口Tool
    """

    def __init__(self, parent):
        super(LdDrawAgent, self).__init__(parent)
        options = [
            {"name": "gro_list", "type": "infile", "format": "dna_evolution.ld_graph"}  # 传入分组路径
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gro_list"):
            raise OptionError("gro_list不存在")

    def set_resource(self):
        self._cpu = 2
        self._memory = "13G"

    def end(self):
        super(LdDrawAgent, self).end()


class LdDrawTool(Tool):
    def __init__(self, config):
        super(LdDrawTool, self).__init__(config)
        # self.R_path = 'program/R-3.3.1/bin/Rscript'
        self.set_environ(
            LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        self.R_path = 'program/miniconda3/envs/R_v4/bin/Rscript'
        self.ld_graph = self.config.PACKAGE_DIR + "/dna_evolution/ld-decay.R"

    def run_ld_graph(self):
        """
        用于LD拟合曲线的绘制

        """
        cmd = "{} {} --infile {} --outfile {}".format(
            self.R_path, self.ld_graph, self.option("gro_list").prop["path"], self.output_dir + '/ld')
        command = self.add_command("ld_graph", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("R脚本ld_decay运行完成")
        else:
            self.set_error("R脚本ld_decay运行失败")

    def run(self):
        super(LdDrawTool, self).run()
        self.run_ld_graph()
        self.end()
