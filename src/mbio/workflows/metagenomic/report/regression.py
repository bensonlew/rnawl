# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class RegressionWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RegressionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_type", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "tax_anno_table", "type": "infile", "format": "meta.profile"},
            {"name": "func_anno_table", "type": "infile", "format": "meta.profile"},
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "func_level", "type": "string", "default": ""},
            {"name": "tax_level", "type": "string", "default": ""},
            {"name": "diversity_type", "type": "string"},
            {"name": "diversity_analysis_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "distance_type", "type": "string"},
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.regression = self.add_tool('statistical.regression')
        self.regression_calculation = self.add_module('statistical.regression_calculation')
        self.sort_func_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_tax_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.get_func_abund_table = self.add_tool('meta.create_abund_table')
        self.get_tax_abund_table = self.add_tool('meta.create_abund_table')
        #self.func_regression_calculation = self.add_module('statistical.regression_calculation')
        self.regression_list = [self.sort_tax_samples, self.sort_func_samples]

    def run_get_tax_abund_table(self):
        self.logger.info("已经运行啦")
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            geneset_table = self.ppm_tool1.option("out_table").prop["path"]
        else:
            geneset_table = self.option('geneset_table')
        self.get_tax_abund_table.set_options({
            'anno_table': self.option('tax_anno_table'),
            'geneset_table': geneset_table,
            'level_type': self.option('tax_level'),
        })
        self.get_tax_abund_table.on("end", self.run_tax_sort_samples)
        self.get_tax_abund_table.run()

    def run_get_func_abund_table(self):
        self.logger.info("已经运行啦")
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            geneset_table = self.ppm_tool2.option("out_table").prop["path"]
        else:
            geneset_table = self.option('geneset_table')
        self.get_func_abund_table.set_options({
            'anno_table': self.option('func_anno_table'),
            'geneset_table': geneset_table,
            'level_type': self.option('func_level'),
        })
        self.get_func_abund_table.on("end", self.run_func_sort_samples)
        self.get_func_abund_table.run()

    def run_tax_sort_samples(self):
        abund_table = self.get_tax_abund_table.option("out_table").prop['path']
        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group'),
        })
        #self.sort_tax_samples.on("end", self.run_tax_regression_calculation)
        self.sort_tax_samples.run()

    def run_func_sort_samples(self):
        if self.option("anno_type") == "go":
            abund_table = self.abundance.option("out_table").prop['path']
        else:
            abund_table = self.get_func_abund_table.option("out_table").prop['path']
        self.sort_func_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group'),
        })
        #self.sort_func_samples.on("end", self.run_func_regression_calculation)
        self.sort_func_samples.run()

    def run_regression_calculation(self):
        tax_abund_table = self.sort_tax_samples.option("out_otu_table").prop['path']
        func_abund_table = self.sort_func_samples.option("out_otu_table").prop['path']
        self.regression_calculation.set_options({
            "tax_abund_file": tax_abund_table,
            "func_abund_file": func_abund_table,
            "diversity_analysis_type": self.option('diversity_analysis_type'),
            "diversity_type": self.option('diversity_type'),
            "distance_type": self.option('distance_type')
        })
        self.regression_calculation.on("end", self.run_regression)
        self.regression_calculation.run()

    def run_regression(self):
        tax_abund_table = self.regression_calculation.option("output_tax").prop['path']
        func_abund_table = self.regression_calculation.option("output_func").prop['path']
        self.regression.set_options({
            "taxon_table": tax_abund_table,
            "func_table": func_abund_table,
        })
        self.regression.on("end", self.set_db)
        self.regression.run()

    def run(self):
        #self.IMPORT_REPORT_DATA = True
        #self.IMPORT_REPORT_DATA_AFTER_END = False
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            self.ppm_tool1 = self.add_tool("meta.cal_ppm")
            self.run_cal_ppm(rely=self.run_get_tax_abund_table, tool_name=self.ppm_tool1, is_run=True)
            if self.option("anno_type") == "go":
                self.run_abundance(self.run_func_sort_samples)
                self.abundance.run()
                self.ppm_tool2 = self.abundance
            else:
                self.ppm_tool2 = self.add_tool("meta.cal_ppm")
                self.run_cal_ppm(rely=self.run_get_func_abund_table, tool_name=self.ppm_tool2, is_run=True)
        else:
            self.run_get_tax_abund_table()
            if self.option("anno_type") == "go":
                self.run_go_abu(rely=self.run_func_sort_samples, is_run=True, level=self.option('func_level'))
            else:
                self.run_get_func_abund_table()
        self.on_rely(self.regression_list, self.run_regression_calculation)
        super(RegressionWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_regression = self.api.api("metagenomic.regression")
        data_path = self.regression.output_dir + "/Regression.data.xls"
        line_path = self.regression.output_dir + "/Regression.message.xls"
        if self.option('diversity_analysis_type') != 'nmds':
            level = 'PC1'
        else:
            level = 'MDS1'
        self.api_regression.add_resgression_site_detail(file_path=data_path, table_id=self.option("main_id"),
                                                        level=level)
        self.api_regression.add_resgression_message_detail(file_path=line_path, table_id=self.option("main_id"),
                                                           level=level)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "regression")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_id"), "regression")
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

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.regression.output_dir))
        result_dir = self.add_upload_dir(self.regression.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "回归分析结果目录", 0, "120208"],
            ["./Regression.data.xls", "xls", "散点图数据", 0, "120209"],
            ["./Regression.message.xls", "xls", "回归曲线信息及R2值", 0, "120210"],
            ["Regression.pdf", "pdf", "物种与功能回归分析图"]
        ])
        super(RegressionWorkflow, self).end()
