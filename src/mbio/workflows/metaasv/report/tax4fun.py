#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author  :  qingchen.zhang


import os
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mbio.packages.metaasv.common_function import link_dir, link_file


class Tax4funWorkflow(Workflow):
    """
    metaasv Tax4Fun 预测分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Tax4funWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的ASV表
            {"name": "input_otu_id", "type": "string"},  # 输入的OTU id
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},  # 输入的group_detail 示例如下
            {"name": "method", "type": "string", "default": ""},
            {"name": "task_id", "type": "string"},
            {"name": "level", "type": "int", "default": 9}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sort_samples = self.add_tool("meta.otu.sort_samples")
        self.tax4fun = self.add_tool("metaasv.tax4fun")
        self.otu_file = self.option("in_otu_table").prop['path']
        self.api_tax4fun = self.api.api("metaasv.tax4fun")
        self.group_table_path = self.option("group").prop['path']

    def run_sort_samples(self):
        """
        排序，按照分组进行合并
        :return:
        """
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.group_table_path,
            "method": self.option("method")
        })
        file_path = os.path.join(self.sort_samples.output_dir, "taxa.table.xls")
        if os.path.isfile(self.otu_file):
            self.otu_file = file_path
        else:
            self.set_error("样本丰度分组计算结果文件不存在%s", variables=(str(file_path)))
        self.sort_samples.run()

    def run_tax4fun(self):
        """
        运行 tax4fun tool
        :return:
        """
        info = self.api_tax4fun.get_database(self.option("task_id"))
        self.logger.info("info: {}".format(info))
        if info['database'] not in ['silva123/16s_bacteria', 'silva123/16s_archaea', 'silva123/16s',
                                    'silva119/16s_bacteria', 'silva119/16s_archaea', 'silva119/16s',
                                    'silva128/16s_archaea', 'silva128/16s_bacteria', 'silva128/16s',
                                    'silva132/16s_archaea', 'silva132/16s_bacteria', 'silva132/16s',
                                    'silva138/16s_archaea', 'silva138/16s_bacteria', 'silva138/16s',
                                    'nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria', 'nt_v20200327/16s',
                                    'greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria',
                                    'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea']:
            raise OptionError("注释数据库%s不是16s细菌或古菌数据库，不能进行Tax4fun分析！", variables=(self.option("database")))
        else:
            self.tax4fun.set_options({
                "in_otu_table": self.otu_file,
                "database": info['database'].encode("utf-8")
            })
            self.tax4fun.run()

    def set_db(self):
        """
        连接结果文件和导入MongoDB
        :return:
        """
        link_dir(self.tax4fun.output_dir, self.output_dir)
        self.logger.info("正在写入mongo数据库")
        files = ["predictions_ko.L1.xls", "predictions_ko.L2.xls", "predictions_ko.L3.xls", "predictions_KO.xls",
                 "predictions_Enzyme.xls"]
        table_name = ["tax4fun_level", "tax4fun_level", "tax4fun_level", "tax4fun_detail", "tax4fun_detail"]
        out_type = ["level1", "level2", "level3", "KO", "enzyme"]
        for i in xrange(len(files)):
            table_path = os.path.join(self.tax4fun.output_dir, files[i])
            if os.path.isfile(table_path):
                self.api_tax4fun.add_tax4fun(table_name[i], table_path, out_type[i], self.option("main_id"))
            else:
                self.set_error("文件不存在%s", variables=(str(table_path)))
        self.end()

    def end(self):
        """
        上传文件和结束
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "KEGG功能预测结果文件",0,""],
            ["predictions_ko.L1.xls", "xls", "KEGG level 1丰度表",0,""],
            ["predictions_ko.L2.xls", "xls", "KEGG level 2丰度表",0,""],
            ["predictions_ko.L3.xls", "xls", "KEGG level 3丰度表",0,""],
            ["predictions_KO.xls", "xls", "KEGG KO丰度表",0,""],
            ["predictions_Enzyme.xls", "xls", "KEGG enzyme丰度表",0,""]
        ])
        super(Tax4funWorkflow, self).end()

    def run(self):
        """
        运行
        :return:
        """
        if self.option("method") == "":
            self.tax4fun.on("end", self.set_db)
            self.run_tax4fun()
        else:
            self.sort_samples.on("end", self.run_tax4fun)
            self.tax4fun.on("end", self.set_db)
            self.run_sort_samples()
        super(Tax4funWorkflow, self).run()
