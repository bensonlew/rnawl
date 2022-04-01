# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __modify__ = '2020/4/02'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir,link_file
import pandas as pd


class AssembleMulAssess2Module(Module):
    """
    多个样品的多个组装结果统计
    """

    def __init__(self, work_id):
        super(AssembleMulAssess2Module, self).__init__(work_id)
        option = [
            {"name": "dir", "type": "string"}
        ]
        self.add_option(option)
        self.run_tools = []
        self.dict = {}

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def run(self):
        super(AssembleMulAssess2Module, self).run()
        self.run_draft_stat()

    def run_draft_stat(self):
        samples = os.listdir(self.option("dir"))
        for sample in samples:
            tool = self.add_module("bacgenome.assemble_mul_assess")
            self.dict[tool] = sample
            scaf_path = self.option("dir") + "/" + sample
            tool.set_options({
                "dir": scaf_path,
                "sample_name": sample,
            })
            self.run_tools.append(tool)
        if len(self.run_tools) >1:
            self.on_rely(self.run_tools, self.set_output)
        elif len(self.run_tools) == 1:
            self.run_tools[0].on("end", self.set_output)
        for tool in self.run_tools:
            tool.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(self.run_tools) >1:
            for tool in self.run_tools:
                link_dir(tool.output_dir, self.output_dir + "/" + self.dict[tool])
        elif len(self.run_tools) == 1:
            link_dir(self.run_tools[0].output_dir, self.output_dir + "/" + self.dict[self.run_tools[0]])
        self.logger.info("设置结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssembleMulAssess2Module, self).end()
