# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

"""组间差异性多组比较检验分析"""

# from biocluster.workflow import Workflow
from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class MultipleModule(Module):
    """
    报告中调用组间差异性分析检验时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MultipleModule, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': " meta.otu.otu_table, toolapps.table,meta.otu.tax_summary_dir"},
            {"name": "group_file", "type": "infile", "format": "toolapps.group_table"},
            # {"name": "type", "type": "string"},
            # {"name": "group_detail", "type": "string"},
            {"name": "test", "type": "string"},
            # {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            # {"name": "params", "type": "string"},
            # {"name": "group_name", "type": "string"},
            {"name": "methor", "type": "string"},
            {"name": "coverage", "type": "float","default":0.95},
            {"name": "ci", "type": "float", "default": 0.05},
            # {"name": "category_name", "type": "string"},
            # {"name": "update_info", "type": "string"},
            # {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.multiple = self.add_tool("statistical.metastat")
        self.sam = self.add_tool("meta.otu.sort_samples_mg")
        self.output_dir = self.multiple.output_dir

    def check_options(self):
        if len(self.option("group_file").prop['group_scheme']) > 1:
            raise OptionError('分组文件不能超过两列')
        if self.option("group_file").group_num(self.option("group_file").prop['group_scheme'][0]) < 3:
            raise OptionError('多组比较组别必须大于等于3')
        group_detail = self.option("group_file").get_group_detail()
        num =[]
        for group_name in group_detail:
            for category in group_detail[group_name]:
                sub_num = len(group_detail[group_name][category])
                num.append(sub_num)
                if sub_num < 2:
                    raise OptionError('每组样品应≥2')
        group_sample = self.option("group_file").prop['sample_name']
        otu_sample = self.option("otu_file").get_sample_info()
        self.logger.info(group_sample)
        self.logger.info(otu_sample)
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

    def run_multiple(self):
        group_name = self.option("group_file").prop['group_scheme'][0]
        if self.option("test") == "anova":
            options = {
                "anova_input": self.sam.option("out_otu_table"),
                "anova_group": self.option("group_file"),
                "anova_correction": self.option("correction"),
                "test": self.option("test"),
                "anova_gname": group_name,
                "anova_methor": self.option("methor"),
                "anova_coverage": self.option("coverage")
            }
        else:
            options = {
                "kru_H_input": self.sam.option("out_otu_table"),
                "kru_H_group": self.option("group_file"),
                "kru_H_correction": self.option("correction"),
                # "kru_H_type": self.option("type"),
                "test": self.option("test"),
                "kru_H_gname": group_name,
                "kru_H_methor": self.option("methor"),
                "kru_H_coverage": self.option("coverage")
            }
        self.multiple.set_options(options)
        self.multiple.on("end", self.end)
        self.multiple.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异多组比较结果目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"],
            [r".*(-).*", "xls", "组间差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值"],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值"]
        ])
        # self.logger.info("------------------")
        # step_main = self.step
        # self.logger.info("spend time is : %s" % step_main.spend_time)
        # (cpu,memory) = self.count_used()
        # self.logger.info("cpu used is : %s" % cpu)
        # self.logger.info("memory used is : %s" % memory)
        # self.logger.info("--------------------")
        super(MultipleModule, self).end()

    def run(self):
        super(MultipleModule, self).run()
        self.logger.info("run diff analysis")
        self.sam.on('end', self.run_multiple)
        self.sort_sample()