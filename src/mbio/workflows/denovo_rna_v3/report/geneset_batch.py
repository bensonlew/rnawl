# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.workflow import Workflow
import os
import time
import pandas as pd
import time


class GenesetBatchWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetBatchWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='genesets', type='string'),
            dict(name='geneset_type', type='string'),
            dict(name="geneset_suffix", type="string", default=None),
            dict(name='diff_id', type='string'),
            dict(name='task_id', type='string'),
            dict(name="update_info", type='string'),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("denovo_rna_v3.export_diff_genesets")

    def run(self):
        self.tool.on("end", self.set_db)
        self.export_diff_genesets()
        super(GenesetBatchWorkflow, self).run()

    def export_diff_genesets(self):
        options = dict(
            diff_id=self.option('diff_id'),
        )
        self.tool.set_options(options)
        self.tool.run()

    def set_db(self):
        geneset_self = self.api.api("denovo_rna_v3.geneset_batch")
        geneset_file = os.path.join(self.tool.output_dir, "genesets")
        geneset_self.add_diff_genesets(self.option('task_id'), geneset_file, self.option('geneset_type'),
                                       self.option('geneset_suffix'), self.option('diff_id'))
        self.end()

    def end(self):
        super(GenesetBatchWorkflow, self).end()
