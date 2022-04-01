# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import pandas as pd
from mbio.packages.bac_comp_genome.common_function import link_dir, link_file


class BivariateBoxWorkflow(Workflow):
    """
    差异倍数
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BivariateBoxWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_table", "type": "infile", "format": "tool_lab.table"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "sample_label", "type": "string", "default": "row"},
            {"name": "main_id", "type": "string"},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.bivariate_box = self.add_tool("tool_lab.bivariate_box")

    def run(self):
        self.run_bivariate_box()
        super(BivariateBoxWorkflow, self).run()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("data_table").is_set:
            raise OptionError("必须输入table文件")
        if not self.option("group_table").is_set:
            raise OptionError("必须输入group文件")

        group = pd.read_table(self.option("group_table").prop["path"], "\t")
        if len(group.columns) != 3:
            raise OptionError("数据表列数必须等于3")
        grouped = group.groupby([group.columns[1], group.columns[2]]).count().reset_index(drop=False)
        for i in grouped[grouped.columns[2]]:
            if i < 4:
                raise OptionError("每个亚组中的样本数必须≥4")

    def run_bivariate_box(self):
        options = {
            "data_table": self.option('data_table').path,
            "sample_label": self.option('sample_label'),
            "group_table": self.option('group_table').path
        }
        self.bivariate_box.set_options(options)
        self.bivariate_box.on("end", self.set_db)
        self.bivariate_box.run()

    def set_db(self):
        self.logger.info("导表开始")
        api_bivariate_box = self.api.api('tool_lab.bivariate_box')
        main_id = self.option("main_id")
        #link_dir(self.bivariate_box.output_dir, self.output_dir)
        box_file = self.bivariate_box.output_dir + "/box_file.txt"
        scatter_file = self.bivariate_box.output_dir + "/scatter_file.txt"
        api_bivariate_box.add_bivariate_box_detail(main_id, box_file=box_file, scatter_file=scatter_file)

        self.end()
        self.logger.info("导表结束")

    def end(self):
        os.link(self.bivariate_box.output_dir + "/box_file.txt", self.output_dir + "/box_data.txt")
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "簇状箱形图结果输出目录"],
            ["./box_data.txt", "txt", "箱形图数据"],

        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BivariateBoxWorkflow, self).end()