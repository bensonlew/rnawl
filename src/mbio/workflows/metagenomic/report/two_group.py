# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

"""组间差异性两组比较检验分析"""

from report import ReportWorkflow
import os
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class TwoGroupWorkflow(CommTableWorkflow):
    """
    报告中调用组间差异性分析检验时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TwoGroupWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "type", "type": "string", "default": "two.side"},
            {"name": "test", "type": "string"},
            # {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "ci", "type": "float", "default": 0.05},
            {"name": "group_name", "type": "string"},
            {"name": "coverage", "type": "float"},
            # {"name": "params", "type": "string"},
            # {"name": "category_name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.two_group = self.add_tool("statistical.metastat")

    def run_two_group(self):
        if not self.option('otu_file').is_set:
            self.option("otu_file").set_path(self.abundance.option('out_table').prop['path'])
        if self.option("test") == "student":
            options = {
                "student_input": self.option("otu_file"),
                "student_group": self.option("group_file"),
                "student_ci": self.option("ci"),
                "student_correction": self.option("correction"),
                "student_type": self.option("type"),
                "test": self.option("test"),
                "student_gname": self.option("group_name"),
                "student_coverage": self.option("coverage")
            }
        elif self.option("test") == "mann":
            options = {
                "mann_input": self.option("otu_file"),
                "mann_ci": self.option("ci"),
                "mann_group": self.option("group_file"),
                "mann_correction": self.option("correction"),
                "mann_type": self.option("type"),
                "test": self.option("test"),
                "mann_gname": self.option("group_name"),
                "mann_coverage": self.option("coverage")
            }
        elif self.option("test") == "signal":
            options = {
                "signal_input": self.option("otu_file"),
                "signal_ci": self.option("ci"),
                "signal_group": self.option("group_file"),
                "signal_correction": self.option("correction"),
                "signal_type": self.option("type"),
                "test": self.option("test"),
                "signal_gname": self.option("group_name"),
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
                "welch_gname": self.option("group_name"),
                "welch_coverage": self.option("coverage")
            }
        self.two_group.set_options(options)
        self.two_group.on("end", self.set_db)
        # self.output_dir = self.two_group.output_dir modify by ghd @ 20181217
        self.two_group.run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "两组比较结果目录", 0, "120185"],
            ["DiffTwoGroup.pdf", "pdf", "两组差异检验柱形图"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "差异显著性比较结果表，包括均值，标准差，p值", 0, "120186"],
            [r".*_CI\.xls", "xls", "差异显著性比较,两组比较置信区间值以及效应量", 0, "120188"],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, "120187"]
        ])
        super(TwoGroupWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        self.output_dir = self.two_group.output_dir  # 防止two_group路径因重运行，导致错误 by ghd @ 20181217
        api_two_group = self.api.api("metagenomic.metastat")
        stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
        ci_path = self.output_dir + '/' + self.option("test") + '_CI.xls'
        bar_path = self.two_group.work_dir + '/' + self.option("test") + '_plot_group_bar.xls'
        if not os.path.isfile(stat_path):
            self.set_error("找不到报告文件:%s", variables=(stat_path), code="12805201")
        if not os.path.isfile(boxfile_path):
            self.set_error("找不到报告文件:%s", variables=(boxfile_path), code="12805202")
        if not os.path.isfile(ci_path):
            self.set_error("找不到报告文件:%s", variables=(ci_path), code="12805203")
        # params = eval(self.option("params"))
        api_two_group.add_metastat_detail(statfile=stat_path, cifiles=[ci_path], check_type='two_group', metastat_id=self.option('main_id'))
        api_two_group.add_metastat_bar(bar_path, self.option('main_id'))
        api_two_group.add_metastat_box(boxfile_path, self.option('main_id'))
        api_two_group.update_metastat(self.option('main_id'), stat_path, ci_path, 'two_group')
        # api_two_group.add_species_difference_check_detail(statfile=stat_path, cifiles=[ci_path], table_id=self.option('main_id'), level=self.option("level"), check_type='two_group', params=self.option("params"), category_name=self.option('category_name'), group_id=params["group_id"], from_otu_table=params["otu_id"], major=False, posthoc=None)
        # api_two_group.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'))
        # print bar_path
        # api_two_group.add_species_difference_check_barplot(bar_path, self.option('main_id'))
        # api_two_group.update_species_difference_check(self.option('main_id'), stat_path, ci_path, 'twogroup')
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
            self.run_two_group()
        else:
            self.logger.info("run abundance")
            #self.run_abundance(self.run_two_group)
            self.run_abundance(self.run_two_group)
            self.abundance.run()
        super(TwoGroupWorkflow, self).run()
