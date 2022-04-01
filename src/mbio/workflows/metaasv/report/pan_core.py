# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import shutil
import os
from biocluster.workflow import Workflow
import re
from mbio.packages.metaasv.common_function import link_dir


class PanCoreWorkflow(Workflow):
    """
    pan_core ASV计算
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PanCoreWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", 'format': "metaasv.otu_table"},#输入ASV表
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},#输入Group表
            {"name": "update_info", "type": "string"},#输入Group表，框架需要
            {"name": "group_detail", "type": "string"},
            {"name": "samples", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.pan_core = self.add_tool("metaasv.pan_core")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.samples = re.split(',', self.option("samples"))

    def check_option(self):
        """
        参数二次检查
        :return:
        """
        self.no_zero_otu = os.path.join(self.work_dir, "input_asv.xls")
        my_sps = self.samples
        self.option("in_otu_table").sub_otu_sample(my_sps, self.no_zero_otu)
        num_lines = sum(1 for line in open(self.no_zero_otu))
        if num_lines < 5:
            self.set_error("ASV表里的ASV数目小于5个！请更换ASV表或者选择更低级别的分类水平！")

    def run_pan_core(self):
        """
        运行tool
        :param no_zero_otu:
        :return:
        """
        self.no_zero_otu = os.path.join(self.work_dir, "input_asv.xls")
        my_sps = self.samples
        self.option("in_otu_table").sub_otu_sample(my_sps, self.no_zero_otu)
        num_lines = sum(1 for line in open(self.no_zero_otu))
        if num_lines < 5:
            self.set_error("ASV表里的ASV数目小于5个！请更换ASV表或者选择更低级别的分类水平！")
        if self.option("group_table").prop["is_empty"]:
            options = {
                "in_otu_table": self.no_zero_otu
            }
        else:
            options = {
                "in_otu_table": self.no_zero_otu,
                "group_table": self.option("group_table")
            }
        self.pan_core.set_options(options)
        self.pan_core.on('end', self.set_db)
        self.pan_core.run()

    def set_db(self):
        """
        导入MongoDB数据
        :return:
        """
        self.logger.info("正在开始链接结果文件目录!")
        link_dir(self.pan_core.output_dir, self.output_dir)
        self.logger.info("正在写入mongo数据库")
        api_pan_core = self.api.api("metaasv.pan_core")
        if self.option("main_id") != "":
            main_id = self.option('main_id')
        else:
            main_id = api_pan_core.add_pan_core(params,group_id)
        pan_path = self.pan_core.option("pan_otu_table").prop['path']
        core_path = self.pan_core.option("core_otu_table").prop['path']
        api_pan_core.add_pan_core_detail(pan_path, main_id, "pan")
        api_pan_core.add_pan_core_detail(core_path, main_id, "core")
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        self.run_pan_core()
        super(PanCoreWorkflow, self).run()

    def end(self):
        """
        结束上传结果文件
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Pan/Core结果目录", 0, ""],
            ["Core.richness.xls", "xls", "Core 表格", 0, ""],
            ["Pan.richness.xls", "xls", "Pan 表格", 0, ""]
        ])
        print self.get_upload_files()
        super(PanCoreWorkflow, self).end()
