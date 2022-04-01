# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""组间差异两样品比较检验分析"""

# from biocluster.workflow import Workflow
from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class TwoSampleModule(Module):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoSampleModule, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "test", "type": "string"},
            # {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "methor", "type": "string"},
            {"name": "coverage", "type": "float"}
        ]
        self.add_option(options)
        self.two_sample = self.add_tool("statistical.metastat")

    def check_options(self):
        otu_sample = self.option("otu_file").get_sample_info()
        if self.option("sample1") not in otu_sample:
            raise OptionError('sample1:{}不存在于丰度文件中'.format(self.option("sample1")))
        if self.option("sample2") not in otu_sample:
            raise OptionError('sample2:{}不存在于丰度文件中'.format(self.option("sample2")))

    def run_two_sample(self):
        if self.option("test") == "chi":
            options = {
                "chi_input": self.option("otu_file"),
                "chi_sample1": self.option("sample1"),
                "chi_sample2": self.option("sample2"),
                "chi_correction": self.option("correction"),
                "test": self.option("test"),
                "chi_methor": self.option("methor"),
                "chi_coverage": self.option("coverage")
            }
        else:
            options = {
                "fisher_input": self.option("otu_file"),
                "fisher_ci": self.option("ci"),
                "fisher_sample1": self.option("sample1"),
                "fisher_sample2": self.option("sample2"),
                "fisher_correction": self.option("correction"),
                "test": self.option("test"),
                "fisher_type": self.option("type"),
                "fisher_methor": self.option("methor"),
                "fisher_coverage": self.option("coverage")
            }
        self.two_sample.set_options(options)
        self.output_dir = self.two_sample.output_dir
        self.two_sample.on("end", self.end)
        self.two_sample.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异两样本比较结果目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"]
        ])
        super(TwoSampleModule, self).end()

    def run(self):
        super(TwoSampleModule, self).run()
        self.logger.info("run diff analysis")
        self.run_two_sample()