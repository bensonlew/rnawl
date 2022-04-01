# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modify by zhujuan 2018.02.02 因mantel检验的环境因子矩阵计算tool换成了distance_calc.py file类型检查跟着修改成"meta.otu.otu_table"

"""报告计算mantel检验"""

import os
from biocluster.workflow import Workflow
from comm_table import CommTableWorkflow


class MantelTestWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MantelTestWorkflow, self).__init__(wsheet_object)
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
            {"name": "env_method", "type": "string"},
            {"name": "env_file", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "units", "type": "string"},  # partial factor
            {"name": "abund_method", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.mantel = self.add_module('statistical.mantel_test')

    def run_sort_samples(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.logger.info(abund_table)
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("group"),
        })
        self.sort_samples.on("end", self.run_mantel_test)
        self.sort_samples.run()

    def run_mantel_test(self):
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        num_lines = open(abund_table, 'r').readlines()
        if len(num_lines) < 3:
            self.set_error('丰度表数据少于2行，请重新设置参数!', code="12801901")
        options = {
            "otutable": self.sort_samples.option("out_otu_table"),
            "factor": self.option("env_file"),
            "otumatrixtype": self.option('abund_method'),
            "factormatrixtype": self.option('env_method')
        }
        if self.option('units'):
            options['partial_factor'] = self.option('units')
        self.mantel.set_options(options)
        self.mantel.on('end', self.set_db)
        self.mantel.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        self.logger.info(self.mantel.output_dir)
        api_mantel_test = self.api.api("metagenomic.mantel_test")
        file_path = ""
        if os.path.exists(self.mantel.output_dir + "/Discompare/mantel_results.txt"):
            file_path = self.mantel.output_dir + "/Discompare/mantel_results.txt"
        if os.path.exists(self.mantel.output_dir + "/Discompare/partial_mantel_results.txt"):
            file_path = self.mantel.output_dir + "/Discompare/partial_mantel_results.txt"
        api_mantel_test.add_mantel_test_detail(file_path, self.option("main_id"))
        self.end()

    def run(self):
        if self.option("abund_file").is_set:
            self.run_sort_samples()
        else:
            #self.run_get_abund_table()
            self.run_abundance(self.run_sort_samples)
            self.abundance.run()
        self.output_dir = self.mantel.output_dir
        super(MantelTestWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "MantelTest分析结果目录", 0, "120194"],
            ["./Discompare", "", "Mantel_Test分析结果目录", 0, "120197"],
            ["./Discompare/mantel_results.txt", "txt", "Mantel_Test分析结果表", 0, "120198"],
            ["./Discompare/partial_mantel_results.txt", "txt", "Partial_Mantel_Test分析结果表", 0, "120199"],
            ["./Facdistance", "", "环境因子矩阵结果目录", 0, "120195"],
            ["./Facdistance/factor_out.xls", "xls", "环境因子矩阵结果表", 0, "120196"],
            ["./Distance", "", "群落/功能矩阵结果目录", 0, "120200"],
            ["./partial", "", "限制环境因子矩阵结果目录", 0, "120202"],
            ["./partial/factor_out.xls", "xls", "限制环境因子矩阵结果表", 0, "120203"]
        ])
        regexps = [
            [r'Distance/%s.*\_taxa.table.xls' % self.option('abund_method'), 'xls', '群落/功能矩阵结果表', 0, "120201"]
        ]
        result_dir.add_regexp_rules(regexps)
        super(MantelTestWorkflow, self).end()
