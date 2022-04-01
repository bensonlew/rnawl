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
from mbio.packages.prok_rna.chart import Chart


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

    def run(self):
        self.start_listener()
        self.fire("start")
        self.set_db()
        self.end()

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因集Venn分析结果目录"],
            ['venn.pdf', "", "venn分析结果图"],
        ])
        super(GenesetVennWorkflow, self).end()

    def set_db(self):
        all_exp = self.api.api("prok_rna.all_exp")
        all_exp.add_venn_tt(
            self.option('main_id'),
        )

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.chart_geneset_venn(self.option("geneset_id"))
        chart.to_pdf()

        # move pdf
        if os.path.exists(os.path.join(self.work_dir, 'geneset.venn.pdf')):
            self.move_pdf(os.path.join(self.work_dir, 'geneset.venn.pdf'), os.path.join(self.output_dir, 'venn.pdf'))
