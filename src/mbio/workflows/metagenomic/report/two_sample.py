# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""组间差异两样品比较检验分析"""

# from biocluster.workflow import Workflow
from report import ReportWorkflow
import os
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class TwoSampleWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoSampleWorkflow, self).__init__(wsheet_object)
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
            {"name": "coverage", "type": "float"},
            # {"name": "params", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_sample = self.add_tool("statistical.metastat")

    def run_two_sample(self):
        if not self.option('otu_file').is_set:
            self.option("otu_file").set_path(self.abundance.option('out_table').prop['path'])
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
        # self.output_dir = self.two_sample.output_dir modify by ghd @ 20181217
        self.two_sample.on("end", self.set_db)
        self.two_sample.run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "两样本比较结果目录", 0, "120189"],
            ["DiffTwoSample.pdf", "pdf", "两样本差异检验柱形图"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "差异显著性比较结果表，包含均值、标准差、p值", 0, "120190"],
            [r".*_CI\.xls", "xls", "差异显著性比较,两样本比较的置信区间值以及效应量", 0, "120191"]
        ])
        super(TwoSampleWorkflow, self).end()

    def set_db(self):
        self.output_dir = self.two_sample.output_dir  # 防止two_group路径因重运行，导致错误 by ghd @ 20181217
        api_two_sample = self.api.api("metagenomic.metastat")
        two_sample_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        ci_path = self.output_dir + '/' + self.option("test") + '_CI.xls'
        if not os.path.isfile(two_sample_path):
            self.set_error("找不到报告文件:%s", variables=(two_sample_path), code="12802701")
        if not os.path.isfile(ci_path):
            self.set_error("找不到报告文件:%s", variables=(ci_path), code="12802702")
        api_two_sample.add_metastat_detail(statfile=two_sample_path, cifiles=[ci_path], check_type='two_sample', metastat_id=self.option('main_id'))
        api_two_sample.update_metastat(self.option('main_id'), two_sample_path, ci_path, 'two_sample')
        # api_two_sample.add_species_difference_check_detail(statfile=two_sample_path, cifiles=[ci_path],
        #                                                    table_id=self.option('main_id'), level=self.option("level"),
        #                                                    check_type='two_sample', params=self.option("params"),
        #                                                    group_id=None, from_otu_table=params["otu_id"], major=False,
        #                                                    posthoc=None)
        # api_two_sample.update_species_difference_check(self.option('main_id'), two_sample_path, ci_path, 'twosample')
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "metastat")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_id"), "metastat")
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        if self.option('otu_file').is_set:
            self.logger.info("run diff analysis")
            self.run_two_sample()
        else:
            self.logger.info("run abundance")
            #self.run_abundance(self.run_two_sample)
            self.run_abundance(self.run_two_sample)
            self.abundance.run()
        super(TwoSampleWorkflow, self).run()
