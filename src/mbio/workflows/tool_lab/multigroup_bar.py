# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import pandas as pd
from mbio.packages.bac_comp_genome.common_function import link_dir, link_file


class MultigroupBarWorkflow(Workflow):
    """
    分组柱形图
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MultigroupBarWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "bar_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "taxonomy_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "combination", "type": "string", "default": "none"},
            {"name": "error_bar", "type": "string", "default": "sd"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.multigroup_bar = self.add_tool("tool_lab.multigroup_bar")

    def run(self):
        self.run_multigroup_bar()
        super(MultigroupBarWorkflow, self).run()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bar_table").is_set:
            raise OptionError("必须输入table文件")
        if not self.option("taxonomy_table").is_set:
            raise OptionError("必须输入taxonomy文件")

        f = pd.read_csv(self.option("bar_table").prop["path"], "\t")
        if self.option("combination") == "average":
            if len(f.columns) < 3:
                raise OptionError("进行组内合并，数据表列数必须大于等于3")
        else:
            if len(f.columns) != 2:
                raise OptionError("不进行组内合并，数据表列数必须等于2")

    def run_multigroup_bar(self):
        options = {
            "bar_table": self.option('bar_table').path,
            "taxonomy_table": self.option('taxonomy_table').path,
            "combination": self.option('combination'),
            "error_bar": self.option('error_bar')
        }
        self.multigroup_bar.set_options(options)
        self.multigroup_bar.on("end", self.set_db)
        self.multigroup_bar.run()

    def set_db(self):
        self.logger.info("导表开始")
        api_multigroup_bar = self.api.api('tool_lab.multigroup_bar')

        result_file = self.multigroup_bar.output_dir + "/result_file.txt"
        combination = self.option('combination')
        api_multigroup_bar.add_multigroup_bar_detail(self.option('main_id'), result_file=result_file, combination=combination)

        self.end()
        self.logger.info("导表结束")

    def end(self):
        if self.option("combination") == "average":

            result_file = self.multigroup_bar.output_dir + "/result_file.txt"
            group_file = pd.read_table(result_file, sep="\t")
            group1 = group_file.ix[:, :-1]
            group1.to_csv(self.output_dir + "/group_result.txt", sep='\t', index=False)
            # os.link(self.multigroup_bar.output_dir + "/group_result.txt", self.output_dir + "/group_result.txt")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "分类柱形图结果输出目录"],
                ["./group_result.txt", "txt", "分类柱形图结果文件"]

            ])
            result_dir.add_regexp_rules([
                ["", "", ""]
            ])
        else:
            pass
        super(MultigroupBarWorkflow, self).end()
