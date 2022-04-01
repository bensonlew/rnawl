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


class O2plsdaWorkflow(Workflow):
    """
    O2plsda分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(O2plsdaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "x_data", "type": "infile", "format": "sequence.profile_table"},
            {"name": "y_data", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group", "type": "infile", "format": "sequence.profile_table"},
            {"name": "oxoy", "type": "string", 'default': "Microbe;Metabolome"},
            {"name": "scale", "type": "string", 'default': "UV"},
            {"name": "x_method", "type": "string", 'default': "log1"},
            {"name": "y_method", "type": "string", 'default': "log1"},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.run_o2plsda()
        super(O2plsdaWorkflow, self).run()

    def rename_name(self,file ,file2, file3, type):
        with open(file, "r") as f, open(file2, "w") as g, open(file3, "w") as s:
            lines = f.readlines()
            num = 0
            s.write(lines[0])
            for line in lines[1:]:
                num += 1
                lin = line.strip().split("\t")
                g.write("{}\t{}\n".format(type + str(num), lin[0]))
                lin[0] = type + str(num)
                s.write("\t".join(lin) + "\n")

    def run_o2plsda(self):
        self.o2plsda = self.add_tool("tool_lab.o2plsda")
        self.rename_name(self.option('x_data').prop['path'], self.work_dir+"/x_data.name", self.work_dir+"/x_data.xls","bsa")
        self.rename_name(self.option('y_data').prop['path'], self.work_dir+"/y_data.name", self.work_dir+"/y_data.xls","tax")
        self.o2plsda.set_options({
            'x_list':self.work_dir+"/x_data.name",
            'y_list':self.work_dir + "/y_data.name",
            'x_data': self.work_dir+"/x_data.xls",
            'y_data': self.work_dir+"/y_data.xls",
            'group': self.option('group'),
            'oxoy': self.option('oxoy'),
            'scale': self.option('scale'),
            "x_method": self.option('x_method'),
            'y_method': self.option('y_method')
        })
        self.o2plsda.on('end', self.set_db)
        self.o2plsda.run()

    def set_db(self):
        """
        导表
        """
        o2plsda = self.api.api('tool_lab.o2plsda')
        path1 =self.o2plsda.output_dir + "/fit_summary.xls"
        path2 = self.o2plsda.output_dir + "/O2PLS_Loadings.xls"
        path3 = self.o2plsda.output_dir + "/O2PLS_Scores.xls"
        o2plsda.add_o2plsda_detail(ObjectId(self.option("main_id")), path1, path2, path3,self.option('oxoy'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.o2plsda.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "O2PLS结果输出目录"],
            ["./fit_summary.xls", "xls", "O2PLS统计情表"],
            ["./O2PLS_Loadings.xls", "xls", "O2PLS分析Loadings信息表"],
            ["./O2PLS_Scores.xls", "xls", "O2PLS分析Scores统计表"],
        ])
        super(O2plsdaWorkflow, self).end()