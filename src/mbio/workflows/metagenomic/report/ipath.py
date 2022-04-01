# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2018/12/12'

from report import ReportWorkflow
import os
import json


class IpathWorkflow(ReportWorkflow):
    """
    宏基因组ipath分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(IpathWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "color_style", "type": "string", "default": "default"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath_tool = self.add_tool('metagenomic.diff.ipath')

    def run(self):
        self.run_abundance(self.run_ipath)
        super(IpathWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_path = self.api.api('metagenomic.common_api')
        self.logger.info("开始进行导表")
        group = json.loads(self.option('group_detail'))
        update_dic = {
            "result_dir": self.sheet.output,
            "table_columns": ','.join(group.keys())
        }
        api_path.add_main_detail(self.ipath_tool.output_dir + "/ipath_input.xls", "ipath_detail",
                                 self.option("main_id"), "ko,color,wide", has_head=False, main_name="ipath_id",
                                 main_table="ipath", update_dic=update_dic)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.ipath_tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "iPath代谢通路分析结果",0,"120285"],
            ["Biosynthesis_of_secondary_metabolities.svg", "svg",
             "Biosynthesis of secondary metabolities"],
            ["Metabolic_pathways.svg", "svg", "Metabolic pathways",0,"120286"],
            ["Regulatory_pathways.svg", "svg", "Regulatory pathways",0,"120287"],
            ["ipath_input.xls", "", "ipath_input file",0,"120288"]
        ])
        super(IpathWorkflow, self).end()

    def run_ipath(self):
        opts = {
            "ko_profile": self.abundance.option('out_table').path,
            "group_table": self.option("group_file"),
            "color_style": self.option("color_style")
        }
        self.ipath_tool.set_options(opts)
        self.ipath_tool.on("end", self.set_db)
        self.ipath_tool.run()
