# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import pandas as pd
from mbio.packages.bac_comp_genome.common_function import link_dir, link_file


class AbundanceBubbleWorkflow(Workflow):
    """
    差异倍数
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AbundanceBubbleWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "taxonomy_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "sample_label", "type": "string", "default": "col"},
            {"name": "combination", "type": "string", "default": "none"},
            {"name": "taxonomy", "type": "string", "default": "none"},
            {"name": "normalize", "type": "string", "default": "none"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.abundance_bubble = self.add_tool("tool_lab.abundance_bubble")

    def run(self):
        self.run_abundance_bubble()
        super(AbundanceBubbleWorkflow, self).run()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("data_table").is_set:
            raise OptionError("必须输入table文件")
        if self.option("combination") != "none":
            if not self.option("group_table").is_set:
                raise OptionError("必须输入group文件")
        if self.option("taxonomy") != "none":
            if not self.option("taxonomy_table").is_set:
                raise OptionError("必须输入taxonomy文件")
        if self.option("normalize") not in ["none", "relative", "log10", "z_score"]:
            raise OptionError("选择是否标准化或标准化方式")

    def run_abundance_bubble(self):
        options = {
            "data_table": self.option('data_table').path,
            "sample_label": self.option('sample_label'),
            "combination": self.option('combination'),
            "group_table": self.option('group_table').path,
            "taxonomy": self.option('taxonomy'),
            "taxonomy_table": self.option('taxonomy_table').path,
            "normalize": self.option('normalize')
        }
        self.abundance_bubble.set_options(options)
        self.abundance_bubble.on("end", self.set_db)
        self.abundance_bubble.run()

    def set_db(self):
        self.logger.info("导表开始")
        abundance_bubble = self.api.api('tool_lab.abundance_bubble')
        abundance_file = self.abundance_bubble.output_dir + "/normalized_result.txt"
        if self.option('taxonomy') != "none":
            taxonomy_file = self.option('taxonomy_table').path
            abundance_bubble.add_abundance_bubble_detail2(self.option("main_id"), abundance_table=abundance_file, taxonomy_table=taxonomy_file)
        else:
            abundance_bubble.add_abundance_bubble_detail(self.option("main_id"), abundance_table=abundance_file)
        self.end()
        self.logger.info("导表结束")

    def end(self):
        if self.option('combination') != "none":
            if self.option('normalize') != "none":
                os.link(self.abundance_bubble.output_dir + "/group_result.txt", self.output_dir + "/group_result.txt")
                os.link(self.abundance_bubble.output_dir + "/normalized_result.txt", self.output_dir + "/normalized_result.txt")
            else:
                os.link(self.abundance_bubble.output_dir + "/group_result.txt", self.output_dir + "/group_result.txt")
        else:
            if self.option('normalize') != "none":
                os.link(self.abundance_bubble.output_dir + "/normalized_result.txt", self.output_dir + "/normalized_result.txt")

        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "丰度气泡图结果输出目录"],
            ["./group_result.txt", "txt", "组内平均值"],
            ["./normalized_result.txt", "txt", "标准化结果"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])

        super(AbundanceBubbleWorkflow, self).end()