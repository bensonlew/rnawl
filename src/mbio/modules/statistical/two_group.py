# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

"""组间差异性两组比较检验分析"""

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class TwoGroupModule(Module):
    """
    报告中调用组间差异性分析检验时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoGroupModule, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "toolapps.group_table"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "test", "type": "string"},
            # {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "coverage", "type": "float", "default": 0.95},
            {"name": "ci_method", "type": "string", "default": "none"}
            # {"name": "params", "type": "string"},
            # {"name": "category_name", "type": "string"},
        ]
        self.add_option(options)
        self.two_group = self.add_tool("statistical.metastat")
        self.sam = self.add_tool("meta.otu.sort_samples_mg")

    def check_options(self):
        self.option("group_file").get_info()
        if len(self.option("group_file").prop['group_scheme']) > 1:
            raise OptionError('分组文件不能超过两列')
        if self.option("group_file").group_num(self.option("group_file").prop['group_scheme'][0]) != 2:
            raise OptionError('两组比较组别必须等于2')
        group_detail = self.option("group_file").get_group_detail()
        num =[]
        for group_name in group_detail:
            for category in group_detail[group_name]:
                sub_num = len(group_detail[group_name][category])
                num.append(sub_num)
                if self.option("test") in ["welch","student"]:
                    if sub_num < 3:
                        raise OptionError('检验方式选择welch和student每组样品应≥3')
                if self.option("test") in ["mann"]:
                    if sub_num < 2:  # 改成和多样性云平台一致，4 -> 2
                        raise OptionError('检验方式选择mann每组样品应≥2')
                if self.option("test") in ["signal"]:
                    if sub_num < 6 or sub_num != num[0]:
                        raise OptionError('检验方式选择signal每组样品应≥6且每组样本数相等')
        group_sample = self.option("group_file").prop['sample_name']
        otu_sample = self.option("otu_file").get_sample_info()
        for sample in group_sample:
            if sample in otu_sample:
                pass
            else:
                raise OptionError('分组文件中样品:{}不存在于丰度文件中'.format(sample))

    def sort_sample(self):
        otutable = self.option("otu_file")
        options = {
            'in_otu_table': otutable,
            'group_table': self.option("group_file")
        }
        self.sam.set_options(options)
        self.sam.run()

    def run_two_group(self):
        group_name = self.option("group_file").prop['group_scheme'][0]
        if self.option("test") == "student":
            options = {
                "student_input": self.sam.option("out_otu_table"),
                "student_group": self.option("group_file"),
                "student_ci": self.option("ci"),
                "student_correction": self.option("correction"),
                "student_type": self.option("type"),
                "test": self.option("test"),
                "student_gname": group_name,
                "student_coverage": self.option("coverage")
            }
        elif self.option("test") == "mann":
            options = {
                "mann_input": self.sam.option("out_otu_table"),
                "mann_ci": self.option("ci"),
                "mann_group": self.option("group_file"),
                "mann_correction": self.option("correction"),
                "mann_type": self.option("type"),
                "test": self.option("test"),
                "mann_gname": group_name,
                "mann_coverage": self.option("coverage")
            }
        elif self.option("test") == "signal":
            options = {
                "signal_input": self.sam.option("out_otu_table"),
                "signal_ci": self.option("ci"),
                "signal_group": self.option("group_file"),
                "signal_correction": self.option("correction"),
                "signal_type": self.option("type"),
                "test": self.option("test"),
                "signal_gname": group_name,
                "signal_coverage": self.option("coverage")
            }
        else:
            options = {
                "welch_input": self.option("otu_file"),
                "welch_ci": self.option("ci"),
                "welch_group": self.option("group_file"),
                "welch_correction": self.option("correction"),
                "welch_type": self.option("type"),
                "test": self.option("test"),
                "welch_gname": group_name,
                "welch_coverage": self.option("coverage")
            }
        self.two_group.set_options(options)
        self.two_group.on("end", self.end)
        self.output_dir = self.two_group.output_dir
        self.two_group.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异两组比较结果目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值"]
        ])
        super(TwoGroupModule, self).end()

    def run(self):
        super(TwoGroupModule, self).run()
        self.logger.info("run diff analysis")
        self.sam.on('end', self.run_two_group)
        self.sort_sample()