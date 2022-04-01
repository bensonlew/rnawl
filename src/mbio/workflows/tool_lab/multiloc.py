# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import json
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from mbio.packages.lnc_rna.copy_file import CopyFile


class MultilocWorkflow(Workflow):
    """
    kegg 编辑
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MultilocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "seq", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {"name": "go", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {'name': 'species', 'type': 'string', 'default': "Animals"},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        print options
        self.add_option(options)
        self.revise_infiles()
        self.multiloc = self.add_tool("tool_lab.multiloc")
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("seq").is_set:
            raise OptionError('请输入ko文件')
        return True

    def run(self):
        self.run_multiloc()
        super(MultilocWorkflow, self).run()

    def run_multiloc(self):
        self.logger.info("开始运行blast注释")
        opts = {
            "fa": self.option("seq").prop['path'],
            "species": self.option("species")
        }

        if self.option("go").is_set:
            opts.update({
                "go": self.option("go").prop['path']
            })
        self.multiloc.set_options(opts)
        self.multiloc.on('end', self.set_output)
        self.multiloc.run()

    def set_db(self):
        # pass
        multiloc_api = self.api.api("tool_lab.multiloc")
        multiloc_api.add_annotation_subloc_detail(self.option("main_id"), os.path.join(self.multiloc.output_dir, "multiloc.xls"))
        multiloc_api.add_annotation_subloc_bar(self.option("main_id"), os.path.join(self.multiloc.output_dir, "multiloc_stat.xls"))


        table_dict = {
            "column": [
                {"field": "accession_id", "title": "Accession", "filter": "false", "sort": "false", "type": "string"},
                {"field": "description", "title": "Description", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc1", "title": "Subcellular Loc1", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc1_prob", "title": "Loc1 Probability", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc2", "title": "Subcellular Loc2", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc2_prob", "title": "Loc2 Probability", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc3", "title": "Subcellular Loc3", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc3_prob", "title": "Loc3 Probability", "filter": "false", "sort": "false", "type": "string"}
            ],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        column_dict = {
            "name": "name",
            "data": "value",
            "condition": {"type": "column"}
        }
        column_data = json.dumps(column_dict, sort_keys=True, separators=(',', ':'))


        multiloc_api.update_db_record('multiloc',
                                      query_dict={"main_id": ObjectId(self.option("main_id"))},
                                      update_dict={'status': 'end',
                                                   'table_data': table_info,
                                                   'column_data': column_data
                                      }
        )


    def set_output(self):
        self.set_db()
        for file in os.listdir(self.multiloc.output_dir):
            CopyFile().linkfile(os.path.join(self.multiloc.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "kegg图片编辑结果",0],
        ])
        super(MultilocWorkflow, self).end()
