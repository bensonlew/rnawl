# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import pandas as pd
import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import matplotlib.pyplot as plt


class MlstWorkflow(Workflow):
    """
    单个基因组的MLST分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MlstWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},# 序列文件
            {"name": "species", "type": "string", },# 物种名称
            {'name': 'task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.run_mlst()
        super(MlstWorkflow, self).run()

    def run_mlst(self):
        self.mlst = self.add_tool("tool_lab.mlst")
        sample = ".".join(os.path.basename(self.option('fasta').prop['path']).split(".")[0:-1])
        self.mlst.set_options({
            'fasta': self.option('fasta'),
            "sample_name": sample,
            'species': self.option('species')
        })
        self.mlst.on('end', self.set_db)
        self.mlst.run()

    def set_db(self):
        """
        导表
        """
        sample = ".".join(os.path.basename(self.option('fasta').prop['path']).split(".")[0:-1])
        mlst = self.api.api('tool_lab.mlst')
        path1 =self.mlst.output_dir + "/" +sample + ".mlst.ST.xls"
        path2 = self.mlst.output_dir + "/" + sample + ".mlst.detail.xls"
        mlst.add_mlst_detail(ObjectId(self.option("main_id")), path1, path2, sample)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.mlst.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Mlst结果输出目录"],
            ["./*.mlst.detail.xls", "xls", "MLST分型比对详情表"],
            ["./*.HitMLST.fasta", "xls", "MLST分型比对等位基因信息"],
            ["./*.mlst.ST.xls", "xls", "MLST分型统计表"],
        ])
        result_dir.add_regexp_rules([
            ["./*.mlst.detail.xls", "xls", "MLST分型比对详情表"],
            ["./*.HitMLST.fasta", "xls", "MLST分型比对等位基因信息"],
            ["./*.mlst.ST.xls", "xls", "MLST分型统计表"],
        ])
        super(MlstWorkflow, self).end()
