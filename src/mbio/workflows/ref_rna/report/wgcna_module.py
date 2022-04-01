# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import re


class WgcnaModuleWorkflow(Workflow):
    """
    666
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(WgcnaModuleWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="main_id", type='string'),
            dict(name="mergeCutHeight", type='string'),
            dict(name="power", type="string"),
            dict(name="networkType", type="string"),
            dict(name="minModuleSize", type="string"),
            dict(name="minKMEtoStay", type="string"),
            dict(name="exp_level", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("rna.wgcna.wgcna_module")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(WgcnaModuleWorkflow, self).run()

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
        dump_tool.update_db_record('sg_wgcna_module', self.option('main_id'), output_dir=output_dir)
        # add result info
        exp_matrix, gene_id2gene_name = self.option('exp_matrix').split(";")
        dump_tool.add_module_detail(
            self.tool.work_dir,
            gene_id2gene_name,
            exp_matrix,
            self.option('main_id')
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "wgcna_module"],
        ])
        super(WgcnaModuleWorkflow, self).end()

    def run_tool(self):
        options = dict(
            datExpr=self.option('exp_matrix').split(";")[0],
            mergeCutHeight=float(self.option('mergeCutHeight')),
            power=int(self.option('power')),
            networkType=self.option('networkType'),
            minModuleSize=int(self.option('minModuleSize')),
            minKMEtoStay=float(self.option('minKMEtoStay')),
        )
        self.tool.set_options(options)
        self.tool.run()
