# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

"""组间差异性多组比较检验分析"""

#from biocluster.workflow import Workflow
from report import ReportWorkflow
import os
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class MultipleWorkflow(CommTableWorkflow):
    """
    报告中调用组间差异性分析检验时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MultipleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            # {"name": "type", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "test", "type": "string"},
            # {"name": "level", "type": "int"},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "params", "type": "string"},
            {"name": "group_name", "type": "string"},
            {"name": "methor", "type": "string"},
            {"name": "coverage", "type": "float"},
            # {"name": "category_name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.multiple = self.add_tool("statistical.metastat")
        # self.output_dir = self.multiple.output_dir modify by ghd @ 20181217

    def run_multiple(self):
        if not self.option('otu_file').is_set:
            self.option("otu_file").set_path(self.abundance.option('out_table').prop['path'])
        if self.option("test") == "anova":
            options = {
                "anova_input": self.option("otu_file"),
                "anova_group": self.option("group_file"),
                "anova_correction": self.option("correction"),
                "test": self.option("test"),
                "anova_gname": self.option("group_name"),
                "anova_methor": self.option("methor"),
                "anova_coverage": self.option("coverage")
            }
        else:
            options = {
                "kru_H_input": self.option("otu_file"),
                "kru_H_group": self.option("group_file"),
                "kru_H_correction": self.option("correction"),
                # "kru_H_type": self.option("type"),
                "test": self.option("test"),
                "kru_H_gname": self.option("group_name"),
                "kru_H_methor": self.option("methor"),
                "kru_H_coverage": self.option("coverage")
            }
        self.multiple.set_options(options)
        self.multiple.on("end", self.set_db)
        self.multiple.run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "多组比较结果目录", 0, "120181"],
            ["DiffMultiple.pdf", "pdf", "多组差异检验柱形图"],
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "差异显著性比较结果表，包括均值，标准差，p值", 0, "120183"],
            # [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"],
            [r".*(-).*", "xls", "差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值", 0, "120182"],
            [r".*_boxfile\.xls", "xls", "差异显著性比较箱线图数据", 0, "120184"]
        ])
        # self.logger.info("------------------")
        # step_main = self.step
        # self.logger.info("spend time is : %s" % step_main.spend_time)
        # (cpu,memory) = self.count_used()
        # self.logger.info("cpu used is : %s" % cpu)
        # self.logger.info("memory used is : %s" % memory)
        # self.logger.info("--------------------")
        super(MultipleWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        self.output_dir = self.multiple.output_dir  # 防止two_group路径因重运行，导致错误 by ghd @ 20181217
        api_multiple = self.api.api("metagenomic.metastat")
        stat_path = self.output_dir + '/' + self.option("test") + '_result.xls'
        boxfile_path = self.output_dir + '/' + self.option("test") + '_boxfile.xls'
        bar_path = self.multiple.work_dir + '/' + self.option("test") + '_plot_group_bar.xls'
        if not os.path.isfile(stat_path):
            self.set_error("找不到报告文件:%s", variables=(stat_path), code="12802001")
        if not os.path.isfile(boxfile_path):
            self.set_error("找不到报告文件:%s", variables=(boxfile_path), code="12802002")
        cifiles = []
        for r, d, f in os.walk(self.output_dir):
            for i in f:
                if self.option("methor") in i:
                    ci_path = r + '/' + i
                    if not os.path.isfile(ci_path):
                        self.set_error("找不到报告文件:%s", variables=(ci_path), code="12802003")
                    cifiles.append(ci_path)
        api_multiple.add_metastat_detail(statfile=stat_path, cifiles=cifiles, check_type='multiple', metastat_id=self.option('main_id'), posthoc=self.option("methor"))
        api_multiple.add_metastat_bar(bar_path, self.option('main_id'))
        api_multiple.add_metastat_box(boxfile_path, self.option('main_id'))
        api_multiple.update_detail(self.option('main_id'))
        # api_multiple.add_species_difference_check_detail(statfile=stat_path, cifiles=cifiles, table_id=self.option('main_id'), level=self.option("level"), check_type='multiple', params=self.option("params"), category_name=self.option('category_name'), group_id=params["group_id"], from_otu_table=params["otu_id"], major=False, posthoc=self.option("methor"))
        # api_multiple.add_species_difference_check_boxplot(boxfile_path, self.option('main_id'))
        # api_multiple.add_species_difference_check_barplot(bar_path, self.option('main_id'))
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

    def report_option(self):
        self.logger.info('option otu_file : ' + self.option('otu_file').prop['path'])
        self.logger.info('option group_file : ' + self.option('group_file').prop['path'])
        self.logger.info('option group_detail : ' + self.option('group_detail'))
        self.logger.info('option update_info : ' + self.option('update_info'))
        self.logger.info('option test : ' + self.option('test'))
        level = self.option('level')
        level = str(level)
        self.logger.info('option level : ' + level)
        self.logger.info('option correction : ' + self.option('correction'))
        self.logger.info('option params : ' + self.option('params'))
        self.logger.info('option group_name : ' + self.option('group_name'))
        self.logger.info('option methor : ' + self.option('methor'))
        coverage = self.option('coverage')
        coverage = str(coverage)
        self.logger.info('option coverage : ' + coverage)
        self.logger.info('option category_name : ' + self.option('category_name'))
        self.logger.info('option main_id : ' + self.option('main_id'))

    def run(self):
        if self.option('otu_file').is_set:
            self.logger.info("run diff analysis")
            self.run_multiple()
        else:
            self.logger.info("run abundance")
            #self.run_abundance(self.run_multiple)
            self.run_abundance(self.run_multiple)
            self.abundance.run()
        # self.report_option()
        super(MultipleWorkflow, self).run()
