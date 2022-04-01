# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mbio.packages.metaasv.common_function import link_dir, link_file
from biocluster.workflow import Workflow
import os


class TwoSampleWorkflow(Workflow):
    """
    metaasv 物种差异分析两样本比较
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoSampleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "test", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "methor", "type": "string"},
            {"name": "coverage", "type": "float"},
            {"name": "params", "type": "string", "default": ""},
            {"name": "update_info", "type": "string", "default": ""},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.two_sample = self.add_tool("statistical.metastat")

    def run_two_sample(self):
        """
        运行两样本比较 metastat
        :return:
        """
        if self.option("test") == "chi":
            options = {
                "chi_input": self.option("otu_file"),
                "chi_sample1": self.option("sample1"),
                "chi_sample2": self.option("sample2"),
                "chi_correction": self.option("correction"),
                "test": self.option("test"),
                "chi_methor": self.option("methor"),
                "chi_coverage": self.option("coverage"),
                "cpu": 16
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
                "fisher_coverage": self.option("coverage"),
                "cpu": 16
            }
        self.two_sample.set_options(options)
        self.two_sample.on("end", self.set_db)
        self.two_sample.run()

    def end(self):
        """
        结束和上传结果文件
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异两样本比较结果目录", 0, ""]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, ""],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量", 0, ""]
        ])
        super(TwoSampleWorkflow, self).end()

    def set_db(self):
        """
        连接结果文件和导入MongoDB
        :return:
        """
        link_dir(self.two_sample.output_dir, self.output_dir)
        api_two_sample = self.api.api("metaasv.stat_test")
        # params = eval(self.option("params"))
        two_sample_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        ci_path = self.output_dir + '/' + self.option("test") + '_CI.xls'
        if not os.path.isfile(two_sample_path):
            self.logger.error("找不到报告文件:{}".format(two_sample_path))
            self.set_error("找不到报告文件")
        if not os.path.isfile(ci_path):
            self.logger.error("找不到报告文件:{}".format(ci_path))
            self.set_error("找不到报告文件")
        api_two_sample.add_species_difference_check_detail(statfile=two_sample_path, cifiles=[ci_path], table_id=self.option('main_id'), level=self.option("level"), major=False, posthoc=None, correlation_key="compare_id", coll_name="two_sample_detail", main_coll="two_sample")
        api_two_sample.update_species_difference_check(self.option('main_id'), two_sample_path, ci_path, 'twosample')
        self.end()

    def run(self):
        self.run_two_sample()
        super(TwoSampleWorkflow, self).run()
