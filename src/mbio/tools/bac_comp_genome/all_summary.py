# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify : 2019.11.19

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.align.blast.xml2table import xml2table
from mbio.packages.bac_comp_genome.common_function import link_file
import re


class AllSummaryAgent(Agent):
    """
    注释汇总文件和可移动元件文件的汇总
    """

    def __init__(self, parent):
        super(AllSummaryAgent, self).__init__(parent)
        options = [
            {"name": "anno_sum", "type": "infile", "format": "sequence.profile_table"}, #注释的总览表
            {"name": "pre_dir", "type": "string",},  # 前噬菌体的目录
            {"name": "isl_dir", "type": "string"},  # 基因组岛的目录
            {"name": "ant_dir", "type": "string"},  # 次级代谢产物的目录
            {"name": "out", "type": "outfile", "format": "sequence.profile_table"}, #比对对齐的序列文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("anno_sum").is_set:
            raise OptionError("必须输入anno_sum文件！")
        return True

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(AllSummaryAgent, self).end()


class AllSummaryTool(Tool):
    def __init__(self, config):
        super(AllSummaryTool, self).__init__(config)
        self._version = "1.0"
        self.python = "/miniconda2/bin/python"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'
        self.out = self.work_dir + "/all.summary.xls"


    def run(self):
        """
        运行
        :return:
        """
        super(AllSummaryTool, self).run()
        self.run_sum()
        self.set_output()
        self.end()

    def run_sum(self):
        cmd = '{} {}all_summary.py -i {} -pre {} -isl {} -ant {} -o {} '.format(self.python, self.package_path, self.option("anno_sum").prop['path'], self.option("pre_dir"), self.option("isl_dir"), self.option("ant_dir"), self.out)
        command = self.add_command("run_sum", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_sum运行完成！")
        else:
            self.set_error("run_sum运行完成运行出错!")


    def set_output(self):
        link_file(self.out, self.output_dir + "/all.summary.xls")
        self.option("out", self.output_dir + "/all.summary.xls")