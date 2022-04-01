# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

import os
import re
import math
import time
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError

class LineChartWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LineChartWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.table"},
            {"name": "ranks", "type": "string", "default": "row"},  # 行、列标签
            {"name": "line_type", "type": "string", "default": ''},  # single, or multi
            {"name": "group_file", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "group_repetition", "type": "string", "default": "true"},  # 组内重复选项
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.line_chart = self.add_tool("tool_lab.line_chart")

    def check_options(self):
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        if not self.option('group_file').is_set:
            raise OptionError("必须设置输入分组文件")
        return True

    def run_tool(self):  # 运行不同参数
        options = {"data_file": self.option('data_file'),
                   "ranks": self.option('ranks'),
                   "line_type":self.option('line_type'),
                   "group_file": self.option('group_file'),
                   "group_repetition": self.option('group_repetition'),
                  }
        self.line_chart.set_options(options)
        self.line_chart.on("end", self.set_output)
        self.line_chart.run()

    def set_output(self):
        output_link = os.path.join(self.output_dir, "group_result.txt")
        if os.path.exists(output_link):
            os.remove(output_link)
        os.link(os.path.join(self.line_chart.output_dir, "group_result.txt"), output_link)
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_line_chart = self.api.api("tool_lab.line_chart")
        main_id = self.option("main_id")
        if main_id is None:
            main_id = api_line_chart.add_line_chart()
        result_table = self.output_dir + "/group_result.txt"
        if self.option('line_type') == "single":
            api_line_chart.add_line_chart_detail(main_id=ObjectId(main_id), input_table=result_table,
                                                    line_type="single", group_repetition=self.option('group_repetition'))
        else:
            api_line_chart.add_line_chart_detail(main_id=ObjectId(main_id), input_table=result_table,
                                                    line_type="multi", group_repetition=self.option('group_repetition'))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_tool()
        super(LineChartWorkflow, self).run()

    def end(self):
        if self.option('group_repetition') == 'true':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "折线图结果输出目录"],
                ["group_result.txt", "txt","折线图结果文件"]
                ])
            result_dir.add_regexp_rules([
                    ["", "", ""]
                ])
        super(LineChartWorkflow, self).end()

