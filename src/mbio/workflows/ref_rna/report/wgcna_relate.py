# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re


class WgcnaRelateWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaRelateWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_eigengenes', type='string'),
            dict(name="main_id", type='string'),
            dict(name="trait_path", type='string'),
            dict(name="trait_type", type="string"),
            dict(name="corr_method", type="string"),
            dict(name="exp_level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("rna.wgcna.wgcna_relate")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(WgcnaRelateWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        dump_tool = self.api.api("ref_rna.wgcna")
        # save workflow output path
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        output_dir = self.workflow_output
        dump_tool.update_db_record('sg_wgcna_relate', self.option('main_id'), output_dir=output_dir)
        # add result info
        seq_annot = self.option('exp_eigengenes').split(";")[2]
        dump_tool.add_relate_detail(
            self.tool.work_dir,
            seq_annot,
            self.option('main_id')
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "wgcna_relate"],
        ])
        super(WgcnaRelateWorkflow, self).end()

    def run_tool(self):
        options = dict(
            datExpr=self.option('exp_eigengenes').split(";")[0],
            MEs=self.option('exp_eigengenes').split(";")[1],
            traits=self.option('trait_path'),
            corType=self.option('corr_method'),
        )
        self.tool.set_options(options)
        self.tool.run()
