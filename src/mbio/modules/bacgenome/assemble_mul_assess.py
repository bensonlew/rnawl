# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/15' #2020.04.02 @gaohao

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir,link_file
import pandas as pd


class AssembleMulAssessModule(Module):
    """
    单个样品的多个组装结果统计
    """

    def __init__(self, work_id):
        super(AssembleMulAssessModule, self).__init__(work_id)
        option = [
            {"name": "dir", "type": "string"},
            {'name': 'sample_name', "type": "string"},  # 样本名
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
        super(AssembleMulAssessModule, self).run()
        self.run_draft_stat()

    def run_draft_stat(self):
        types = os.listdir(self.option("dir"))
        for i in types:
            tool = self.add_module("bacgenome.assemble_assess")
            self.dict[tool] = i
            scaf_path = self.option("dir") + "/" + i + "/" + self.option("sample_name") + ".scaffold.fna"
            tool.set_options({
                "seq_scaf": scaf_path,
                "sample_name": self.option("sample_name"),
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
                link_dir(tool.output_dir + '/assembly', self.output_dir + "/" + self.dict[tool])
        elif len(self.run_tools) == 1:
            link_dir(self.run_tools[0].output_dir + '/assembly', self.output_dir + "/" + self.dict[self.run_tools[0]])
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
        super(AssembleMulAssessModule, self).end()
