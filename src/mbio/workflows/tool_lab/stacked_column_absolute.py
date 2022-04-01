# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

import os
import re
import math
import time
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId

class StackedColumnAbsoluteWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(StackedColumnAbsoluteWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.table"},
            {"name": "group_merge", "type": "string", "default": "average"},  # 组内合并，average or none
            {"name": "group_file", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.stacked_column_absolute = self.add_tool("tool_lab.stacked_column_absolute")

    def check_options(self):
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        if self.option("group_merge") == "average":# 选择组内合并时，输入分组文件
            if not self.option("group_file").is_set:
                raise OptionError("必须设置输入分组文件")
        return True

    def run_stacked_column(self):
        options = {"data_file": self.option('data_file'),
                   "group_merge":self.option('group_merge')}
        if self.option('group_merge') == "average":
            options.update({"group_file": self.option('group_file')})
        self.stacked_column_absolute.set_options(options)
        self.stacked_column_absolute.on("end", self.set_output)
        self.stacked_column_absolute.run()

    def set_output(self):
        output_link = os.path.join(self.output_dir, "group_result.txt")
        if os.path.exists(output_link):
            os.remove(output_link)
        output_dir = os.path.join(self.stacked_column_absolute.output_dir, "group_result.txt")
        os.link(output_dir, output_link)
        self.set_db()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        self.logger.info("开始导表")
        api_stacked_column_absolute = self.api.api("tool_lab.stacked_column_absolute")
        table = self.output_dir + "/group_result.txt"
        main_id = self.option("main_id")
        if main_id is None:
            main_id = api_stacked_column_absolute.add_stacked_column_absolute()
        api_stacked_column_absolute.add_stacked_column_absolute_detail(result_file=table,main_id=ObjectId(main_id))
        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_stacked_column()
        super(StackedColumnAbsoluteWorkflow, self).run()

    def end(self):
        if self.option('group_merge') == "average":
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./group_result.txt", "txt", "堆叠柱形图结果文件"]
            ])
            result_dir.add_regexp_rules([
                ["", "", ""]
            ])
        super(StackedColumnAbsoluteWorkflow, self).end()

