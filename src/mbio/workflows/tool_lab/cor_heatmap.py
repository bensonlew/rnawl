# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import shutil
from bson.objectid import ObjectId


class CorHeatmapWorkflow(Workflow):
    """
    checkM

    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CorHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table1", "type": "infile", "format": "tool_lab.table"},
            {"name": "table2", "type": "infile", "format": "tool_lab.table"},
            {"name": "method", "type": "string", "default": "pearsonr"},
            {"name": "row_cluster", "type": "string", "default": "none"},
            {"name": "column_cluster", "type": "string", "default": "none"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_corheatmap()
        super(CorHeatmapWorkflow, self).run()

    def run_corheatmap(self):
        self.corheatmap = self.add_tool("tool_lab.cor_heatmap")

        if self.option('row_cluster') == "none":
            row_cluster = ""
        else:
            row_cluster = self.option('row_cluster')
        if self.option('column_cluster') == "none":
            column_cluster = ""
        else:
            column_cluster = self.option('column_cluster')
        options = {
            "table1": self.option('table1'),
            "table2": self.option('table2'),
            "method": self.option('method'),
            "row_cluster": row_cluster,
            "column_cluster": column_cluster,
        }
        self.corheatmap.set_options(options)
        self.corheatmap.on('end', self.set_output)
        self.corheatmap.run()

    def set_output(self):
        for file in os.listdir(self.corheatmap.output_dir):
            os.link(os.path.join(self.corheatmap.output_dir, file), os.path.join(self.output_dir, file))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_corheatmap = self.api.api("tool_lab.cor_heatmap")
        api_corheatmap.add_corheatmap(self.corheatmap.output_dir, main_id=self.option('main_id'))
        self.logger.info("导表结束")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.corheatmap.output_dir)
        super(CorHeatmapWorkflow, self).end()