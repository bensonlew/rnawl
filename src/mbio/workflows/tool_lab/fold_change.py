# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
from mbio.packages.bac_comp_genome.common_function import link_dir, link_file
import pandas as pd


class FoldChangeWorkflow(Workflow):
    """
    差异倍数
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FoldChangeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bar_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "sample_label", "type": "string", "default": "col"},
            {"name": "control_group", "type": "string"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.fold_change = self.add_tool("tool_lab.fold_change")

    def run(self):
        self.run_fold_change()
        super(FoldChangeWorkflow, self).run()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bar_table").is_set:
            raise OptionError("必须输入table文件")
        if not self.option("group_table").is_set:
            raise OptionError("必须输入group文件")
        if not self.option("control_group"):
            raise OptionError("必须输入对照组")

        group = self.option("group_table").prop['path']
        with open(group) as f:
            lines = f.readlines()
            group_list = []
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[1] not in group_list:
                    group_list.append(item[1])
        if len(group_list) != 2:
            raise OptionError("分组数不等于2，无法进行分析")
        if self.option("control_group") not in group_list:
            raise OptionError("分析选择的对照组不存在，请检查后重新运行")

    def run_fold_change(self):
        options = {
            "bar_table": self.option('bar_table').path,
            "group_table": self.option('group_table').path,
            "sample_label": self.option('sample_label'),
            "control_group": self.option('control_group')
        }
        self.fold_change.set_options(options)
        self.fold_change.on("end", self.set_db)
        self.fold_change.run()

    def set_db(self):
        self.logger.info("导表开始")
        api_fold_change = self.api.api('tool_lab.fold_change')
        link_dir(self.fold_change.output_dir, self.output_dir)
        result_file = self.output_dir + "/result.txt"
        api_fold_change.add_fold_change_detail(self.option('main_id'), result_file=result_file)
        self.end()
        self.logger.info("导表结束")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "差异倍数结果输出目录"],
            ["./result.txt", "txt", "差异倍数"]

        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(FoldChangeWorkflow, self).end()