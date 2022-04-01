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



class KeggEditorWorkflow(Workflow):
    """
    kegg 编辑
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(KeggEditorWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input1", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {"name": "input2", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {"name": "color_bg1", "type": "string", "default": "#FF0000"},
            {"name": "color_bg2", "type": "string", "default": "#00FF00"},
            {"name": "color_fg11", "type": "string", "default": "#FF0000"},
            {"name": "color_fg12", "type": "string", "default": "#FF0000"},
            {"name": "color_fg21", "type": "string", "default": "#00FF00"},
            {"name": "color_fg22", "type": "string", "default": "#00FF00"},
            {"name": "table_num", "type": "int", "default": "1"},

            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        print options
        self.add_option(options)
        self.revise_infiles()
        self.kegg_editor = self.add_tool("tool_lab.kegg_editor")
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("input1").is_set:
            raise OptionError('请输入ko文件')
        return True

    def run(self):
        self.run_kegg_editor()
        super(KeggEditorWorkflow, self).run()

    def run_kegg_editor(self):
        self.logger.info("开始运行blast注释")
        opts = {
            "input1": self.option("input1"),
            "input2": self.option("input2"),
            "color_bg1": self.option("color_bg1"),
            "color_bg2": self.option("color_bg2"),
            "color_fg11": self.option("color_fg11"),
            "color_fg12": self.option("color_fg12"),
            "color_fg21": self.option("color_fg21"),
            "color_fg22": self.option("color_fg22")
        }

        self.kegg_editor.set_options(opts)
        self.kegg_editor.on('end', self.set_output)
        self.kegg_editor.run()

    def set_output(self):
        for file in os.listdir(self.kegg_editor.output_dir):
            os.link(os.path.join(self.kegg_editor.output_dir, file), os.path.join(self.output_dir, file))
        self.end()
    def set_db(self):
        pass
        # orf_predict_api = self.api.api("tool_lab.orf_predict")
        # orf_predict_api.add_orf_detail(self.option("main_id"), os.path.join(self.cds_predict.output_dir, "all_predicted.xls"))
        # table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        # orf_predict_api.update_db_record('orf_predict',
        #                                  query_dict={"main_id": ObjectId(self.option("main_id"))},
        #                                  update_dict={'status': 'end',
        #                                                                                      'table_data': table_info})


    def set_output(self):
        self.set_db()
        for file in os.listdir(self.kegg_editor.output_dir):
            CopyFile().linkfile(os.path.join(self.kegg_editor.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "kegg图片编辑结果",0],
        ])
        super(KeggEditorWorkflow, self).end()
