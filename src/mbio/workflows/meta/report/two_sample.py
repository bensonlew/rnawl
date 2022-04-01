# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""组间差异两样品比较检验分析"""

from biocluster.workflow import Workflow
import os
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class TwoSampleWorkflow(Workflow):
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
            {"name": "params", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_sample = self.add_tool("statistical.metastat")

    def run_two_sample(self):
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
        self.output_dir = self.two_sample.output_dir
        self.two_sample.on("end", self.set_db)
        self.two_sample.run()

    def end(self):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种差异两样本比较结果目录", 0, "110143"],
            [r"./两样本比较差异检验柱形图.pdf", "pdf", "两样本比较中具有显著差异的丰度排行前15的物种差异检验柱形图", 0, ""]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值", 0, "110144"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量", 0, "110145"]
        ])
        super(TwoSampleWorkflow, self).end()

    def set_db(self):
        api_two_sample = self.api.stat_test
        params = eval(self.option("params"))
        two_sample_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        ci_path = self.output_dir + '/' + self.option("test") + '_CI.xls'
        if not os.path.isfile(two_sample_path):
            self.logger.error("找不到报告文件:{}".format(two_sample_path))
            self.set_error("找不到报告文件", code="12704201")
        if not os.path.isfile(ci_path):
            self.logger.error("找不到报告文件:{}".format(ci_path))
            self.set_error("找不到报告文件", code="12704201")
        api_two_sample.add_species_difference_check_detail(statfile=two_sample_path, cifiles=[ci_path], table_id=self.option('main_id'), level=self.option("level"), check_type='two_sample', params=self.option("params"), group_id=None, from_otu_table=params["otu_id"], major=False, posthoc=None)
        api_two_sample.update_species_difference_check(self.option('main_id'), two_sample_path, ci_path, 'twosample')
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_species_difference_check")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "species_difference_two_specimen",
                "interaction": 1,
                "main_table": "sg_species_difference_check",
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_two_sample()
        super(TwoSampleWorkflow, self).run()
