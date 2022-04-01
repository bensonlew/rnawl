# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re,os
import time
from biocluster.workflow import Workflow
from mbio.packages.whole_transcriptome.chart.chart import Chart
import glob
from mbio.packages.project_demo.interaction_rerun.interaction_delete import linkfile

class GenesetVennWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/10 GeneSet_venn')

    def run(self):
        self.start_listener()
        self.fire("start")
        self.set_db()
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        chart.chart_geneset_venn(self.option("geneset_id"))
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            linkfile(p, self.output_dir + "/" + os.path.basename(p))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', '基因集Venn分析', 0],
            ['*venn.*pdf', 'pdf', '基因集Venn图', 0],
            ['*upset*.pdf', 'pdf', '基因集Upset图', 0],
        ])
        super(GenesetVennWorkflow, self).end()

    def set_db(self):
        all_exp = self.api.api("whole_transcriptome.all_exp")
        all_exp.add_venn_tt(
            self.option('main_id'),
        )