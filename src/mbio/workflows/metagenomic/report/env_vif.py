# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

"""报告计算VIF方差膨胀因子"""

import os
import re
import glob
from biocluster.workflow import Workflow
from comm_table import CommTableWorkflow
from mbio.packages.meta.common_function import envname_restore


class EnvVifWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EnvVifWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_id", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "viflim", "type": "int", "default": 10},
            {"name": "env_file", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "abund_method", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.vif = self.add_tool('statistical.env_vif')
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")

    def run_sort_samples(self):
        self.logger.info("正常运行啦")
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            self.logger.info("对啦对啦")
            abund_table = self.abundance.option("out_table").prop['path']
        self.logger.info(abund_table)
        # self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("group"),
            "sample_del_warn": "T"
        })
        self.sort_samples.on("end", self.run_env_vif)
        self.sort_samples.run()

    def run_env_vif(self):
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        num_lines = open(abund_table, 'r').readlines()
        if len(num_lines) < 3:
            raise OptionError('丰度表数据少于2行，请重新设置参数!', code="12801501")
        options = {
            "abundtable": self.sort_samples.option("out_otu_table"),
            "envtable": self.option('env_file'),
            "viflim": self.option("viflim")
        }
        self.vif.set_options(options)
        self.vif.on('end', self.set_db)
        self.vif.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        api_env_vif = self.api.api("metagenomic.env_vif")
        env_result = glob.glob(self.vif.output_dir + "/*")
        for file in env_result:
            if re.search("final_.*_vif.txt", file):
                final_path = file
                api_env_vif.add_env_vif_detail(final_path, "after_vif", self.option("main_id"))
            elif re.search("raw_.*_vif.txt", file):
                raw_path = file
                api_env_vif.add_env_vif_detail(raw_path, "before_vif", self.option("main_id"))
        self.end()

    def run(self):
        if self.option("abund_file").is_set:
            self.run_sort_samples()
        else:
            #self.run_get_abund_table()
            self.run_abundance(self.run_sort_samples)
            self.abundance.run()
        self.output_dir = self.vif.output_dir
        super(EnvVifWorkflow, self).run()

    @envname_restore
    def end(self):
        env_result = glob.glob(self.vif.output_dir + "/*")
        result_dir = self.add_upload_dir(self.output_dir)
        if re.search("cca_vif.txt", env_result[0]):
            result_dir.add_relpath_rules([
                [".", "", "VIF方差膨胀因子分析结果目录", 0, "120204"],
                ["./raw_cca_vif.txt", "txt", "筛选前VIF方差膨胀因子分析结果表", 0, "120205"],
                ["./final_cca_vif.txt", "txt", "筛选后VIF方差膨胀因子分析结果表", 0, "120206"],
                ["./DCA.txt", "txt", "判断计算VIF方差膨胀因子分析的方法文件", 0, "120207"]
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "VIF方差膨胀因子分析结果目录", 0, "120204"],
                ["./raw_rda_vif.txt", "txt", "筛选前VIF方差膨胀因子分析结果表", 0, "120205"],
                ["./final_rda_vif.txt", "txt", "筛选后VIF方差膨胀因子分析结果表", 0, "120206"],
                ["./DCA.txt", "txt", "判断计算VIF方差膨胀因子分析的方法文件", 0, "120207"]
            ])
        super(EnvVifWorkflow, self).end()
