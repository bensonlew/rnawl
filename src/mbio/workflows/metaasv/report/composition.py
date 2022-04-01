# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os,re
from biocluster.workflow import Workflow
from mbio.packages.metaasv.common_function import link_dir


class CompositionWorkflow(Workflow):
    """
    metaasv 群落组成分析模块
    barpie、circos图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompositionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "asv_id", "type": "string"},  # 输入的ASV id
            {"name": "level", "type": "string", "default": "9"},  # 输入的OTU level
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},##主表id
            {"name": "group_detail", "type": "string"}, # 输入的group_detail 示例如下
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},##group 表
            {"name": "method", "type": "string", "default": ""},#聚类方法
            {"name": "combine_value", "type": "string", "default": ""},#others合并
            {"name": "graphic_type", "type": "string", "default": ""}  # barpie,circos
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sort_samples = self.add_tool("metaasv.sort_samples")
        self.group_table_path = self.option("group")

    def run_sort_samples(self):
        """
        运行tool-sort_samples_mg
        取全表然后top
        :return:
        """
        self.sort_samples.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "group_table": self.group_table_path,
            "method": self.option("method"),
            "others": float(self.option("combine_value")),
            "fill_zero": "false",
            "abundance_method": "all"
        })
        self.sort_samples.on("end", self.set_db)
        self.sort_samples.run()

    def set_db(self):
        """
        导入MongoDB
        :return:
        """
        link_dir(self.sort_samples.output_dir, self.output_dir)
        self.logger.info("正在写入mongo数据库")
        type = self.option("graphic_type")
        if type in ['barpie']:
            api_barpie = self.api.api("metaasv.barpie")
            composition_psth = os.path.join(self.output_dir, "taxa.percents.table.xls")
            api_barpie.add_sg_otu_detail(composition_psth,
                                         self.option("main_id"))
        elif type in ['circos']:
            api_barpie = self.api.api("metaasv.circos")
            composition_psth = os.path.join(self.output_dir, "taxa.percents.table.xls")
            api_barpie.add_sg_otu_detail(composition_psth,
                                         self.option("main_id"))
        self.end()

    def end(self):
        """
        结束，上传文件
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "物种组成分析结果目录", 0, ""],
            ["taxa.table.xls", "xls", "各样本物种丰度结果表", 0, ""],
            ["taxa.precents.table.xls", "xls", "各样本物种相对丰度结果表", 0, ""]
        ])
        super(CompositionWorkflow, self).end()

    def run(self):
        """
        运行
        :return:
        """
        self.run_sort_samples()
        super(CompositionWorkflow, self).run()
