# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/9'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd


class HiseqSeqStatModule(Module):
    """
    分别对二代数据原始文件或质控后文件做两种碱基统计
    """

    def __init__(self, work_id):
        super(HiseqSeqStatModule, self).__init__(work_id)
        option = [
            {"name": "list", "type": "infile", "format": "bacgenome.simple_file"},
            {"name": "fastq_dir", "type": "infile", "format": "bacgenome.simple_dir"}
        ]
        self.add_option(option)
        self.stat_tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(HiseqSeqStatModule, self).run()
        self.get_list()
        self.run_fastq_stat()
        self.run_fastx()
        self.logger.info(self.stat_tools)
        self.on_rely(self.stat_tools, self.set_output)
        self.logger.info(self.stat_tools)
        for tool in self.stat_tools:
            tool.run()

    def get_list(self):
        data = pd.read_table(self.option("list").prop["path"], header=None, index_col=[1,2])
        data[0] = self.option("fastq_dir").prop["path"] + "/" + data[0]
        data = data.unstack()
        data.columns = data.columns.get_level_values(1)
        self.samples = data.to_dict(orient="index")
        #self.logger.info(self.samples)
        #self.set_error("terminate")

    def run_fastq_stat(self):
        for sample in self.samples:
            fastq_stat = self.add_tool("bacgenome.bac_read_statistic2")
            options = {
                "q1": self.samples[sample]["l"],
                "q2": self.samples[sample]["r"],
                "prefix": sample
            }
            fastq_stat.set_options(options)
            self.stat_tools.append(fastq_stat)

    def run_fastx(self):
        for sample in self.samples:
            for fastq in self.samples[sample].values():
                fastx = self.add_tool("bacgenome.fastx_v2")
                options = {
                    "fastq": fastq
                }
                fastx.set_options(options)
                self.stat_tools.append(fastx)

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.logger.info(len(self.stat_tools))
        # edit output
        self.logger.info("设置注释结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(HiseqSeqStatModule, self).end()
