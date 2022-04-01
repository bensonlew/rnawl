# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modifies 20180409

"""circos 图"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.config import Config


class CgviewWorkflow(Workflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CgviewWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "xml_file", "type": "infile", "format": "bacgenome.cgview_xml"},
            {"name": "task_id", "type": "string"},
            {"name": "specimen_id", "type": "string"},
            {"name": "genome_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "species_name", "type": "string"},
            {"name": "seq_type", "type": "string", "default": "Circular"},  # 开闭环
            {"name": "params", "type": "string"},
            {"name": "main_table_data", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.cgview = self.add_tool("bacgenome.cgview")
        self.file_path = self._sheet.output

    def run_cgview(self):
        if self.option('genome_type') in ['scaffold','Scaffold']:
            options = {
                "sample_name":self.option('specimen_id'),
                "xml_file":self.option('xml_file'),
                "species_name":self.option('species_name'),
                "seq_type_circle": self.option('seq_type')
            }
        else:
            options = {
                "sample_name": self.option('specimen_id'),
                "xml_file": self.option('xml_file'),
                "species_name": self.option('species_name'),
                "seq_type": self.option('genome_type'),
                "seq_type_circle": self.option('seq_type')
            }
        self.cgview.set_options(options)
        self.cgview.on('end', self.set_db)
        self.cgview.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_cgview()
        super(CgviewWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        result = self.cgview.output_dir
        api_cgview = self.api.api('bacgenome.cgview')
        api_cgview.add_cgview(result, main=False, main_id=self.option('main_id'),
                              task_id=self.option('task_id'), update=False, link_path=self.file_path)
        self.end()

    def end(self):
        repaths = [
            [".", "", ""],
        ]
        regexps = [
        ]
        sdir = self.add_upload_dir(self.cgview.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CgviewWorkflow, self).end()
