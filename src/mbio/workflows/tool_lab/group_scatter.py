# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

import os
import time
import pandas as pd
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId


class GroupScatterWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GroupScatterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_file", "type": "infile", "format": "toolapps.group_table"},
            {"name": "ranks", "type": "string", "default": "row"},         # 选择行、列标签
            {"name": "main_id", "type": "string"},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.group_scatter = self.add_tool("tool_lab.group_scatter")

    def check_options(self):
        self.logger.info(self.option('data_file'))
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        if not self.option("group_file").is_set:
            raise OptionError("必须设置输入分组文件")
        return True

    def run_group_scatter(self):
        option = {
                "data_file": self.option('data_file'),
                "group_file": self.option('group_file'),
                "ranks":self.option('ranks')
            }
        self.group_scatter.set_options(option)
        self.group_scatter.on("end", self.set_output)
        self.group_scatter.run()

    def set_output(self):
        output_link = os.path.join(self.output_dir, "group_result.txt")
        if os.path.exists(output_link):
            os.remove(output_link)
        output_dir = os.path.join(self.group_scatter.output_dir, "group_result.txt")
        os.link(output_dir, output_link)
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_group_scatter = self.api.api("tool_lab.group_scatter")
        main_id = self.option("main_id")
        if main_id == None:
            main_id = api_group_scatter.add_group_scatter()
        result_table = self.output_dir + "/group_result.txt"
        api_group_scatter.add_group_scatter_box_detail(result_file=result_table, main_id=ObjectId(main_id))
        group_table = self.option('group_file').path
        if self.option('ranks') == "row":
            data_table = self.option('data_file').path
        else:
            data_table = self.group_scatter.output_dir + "/data_file_t.txt"
        api_group_scatter.add_group_scatter_detail(data_file=data_table, group_file=group_table, main_id=ObjectId(main_id))
        self.end()

    def run(self):
        self.run_group_scatter()
        super(GroupScatterWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./group_result.txt", "txt", "分组散点图结果文件"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(GroupScatterWorkflow, self).end()
