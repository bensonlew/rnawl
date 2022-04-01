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
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
import glob
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


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
                                                        'interaction_results/03 GeneSet/10 Geneset_Venn')

    def run(self):
        self.start_listener()
        self.fire("start")
        time.sleep(1)
        self.get_run_log()
        self.set_db()
        self.end()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_geneset_venn", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        chart.chart_geneset_venn(self.option("geneset_id"))
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "目标基因集Venn分析结果文件", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*venn.pdf', 'txt', '目标基因集Venn图', 0],
            ['*upset.pdf', 'txt', '目标基因集Upset图', 0],
        ])
        super(GenesetVennWorkflow, self).end()

    def set_db(self):
        all_exp = self.api.api("ref_rna_v2.all_exp")
        all_exp.add_venn_tt(
            self.option('main_id'),
        )
