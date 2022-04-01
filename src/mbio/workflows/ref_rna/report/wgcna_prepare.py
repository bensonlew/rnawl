# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re


class WgcnaPrepareWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaPrepareWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="main_id", type='string'),
            dict(name="me", type='string'),
            dict(name="cv", type="string"),
            dict(name="group_dict", type="string"),
            dict(name="group_id", type="string"),
            dict(name="exp_level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("rna.wgcna.wgcna_prepare")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(WgcnaPrepareWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("ref_rna.wgcna")
        # add result info
        dump_tool.add_prepare_detail(self.tool.output_dir, main_id=self.option('main_id'), )
        with open(self.work_dir + '/group_info.txt') as f:
            _ = f.readline()
            group_info = dict()
            for line in f:
                s, g = line.strip().split('\t')
                group_info[s] = g
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        exp_matrix =self.workflow_output + '/exp_matrix_after_filtering.txt'
        dump_tool.update_db_record('sg_wgcna_prepare', self.option('main_id'), group_info=group_info, exp_matrix=exp_matrix)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "wgcna_preprocessing"],
        ])
        super(WgcnaPrepareWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
            me=self.option('me'),
            cv=self.option('cv'),
        )
        self.tool.set_options(options)
        self.tool.run()
