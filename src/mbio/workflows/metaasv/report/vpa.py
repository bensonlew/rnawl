# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir


class VpaWorkflow(Workflow):
    """
    metaasv VPA分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VpaWorkflow, self).__init__(wsheet_object)
        options = [
             {"name": "asv_id", "type": "string"},
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},##ASV表
            {"name": "level", "type": "int", "default": 9},
            {"name":"group","type":"infile","format": "meta.otu.group_table"},##group表
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name":"update_info","type":"string"},
            {"name":"main_id","type":"string"},##主表ID
            {"name": "env_file", "type": "infile", "format": "meta.otu.otu_table"} , #环境因子表
            {"name": "env_group", "type": "infile", "format": "meta.otu.group_table"},##环境因子分组表
            {"name": "env_detail", "type": "string"},
            {"name": "env_labs", "type": "string"}##选择的环境因子
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sort_func_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.vpa = self.add_tool('statistical.vpa')

    def run_func_sort_samples(self):
        """
        排序
        :return:
        """
        abund_table = self.option('otu_table')
        self.sort_func_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group'),
        })
        self.sort_func_samples.on("end", self.run_vpa)
        self.sort_func_samples.run()

    def run_vpa(self):
        """
        VPA分析
        :return:
        """
        otu_table = self.sort_func_samples.option("out_otu_table")
        self.vpa.set_options({
            "species_table":  otu_table,
            "env_table": self.option("env_file"),
            "group_table": self.option('env_group')
        })
        self.vpa.on("end", self.set_db)
        self.vpa.run()


    def run(self):
        """
        运行
        :return:
        """
        self.run_func_sort_samples()
        super(VpaWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_vpa = self.api.api("metaasv.vpa")
        link_dir(self.vpa.output_dir, self.output_dir)
        data_path = self.vpa.output_dir + '/env.R2adj.xls'
        graph_path = self.vpa.work_dir + '/env.plot.xls'
        self.api_vpa.add_vpa_detail(data_path,self.option("main_id"))
        self.api_vpa.add_vpa_graph(graph_path,self.option("main_id"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.vpa.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "VPA结果文件目录", 0, ""],
            ["./env.R2adj.xls", "xls", "R2adj数据", 0, ""],
        ])
        super(VpaWorkflow, self).end()
