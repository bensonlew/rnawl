# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
# from mainapp.controllers.project.meta_controller import MetaController
# from mainapp.libs.param_pack import group_detail_sort
# from bson import ObjectId
# import datetime
# import shutil
import re,os
# import time
import glob
from biocluster.workflow import Workflow
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset

class DiffGenesetVennWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffGenesetVennWorkflow, self).__init__(wsheet_object)
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
                                                        'interaction_results/01 Diff_Express/05DiffExp_Venn')

    def run(self):
        self.start_listener()
        self.fire("start")
        self.set_db()
        self.end()

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["01 Diff_Express", "", "差异基因数据挖掘结果目录", 0],
            ["01 Diff_Express/05DiffExp_Venn", "", "差异基因集venn分析", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "差异基因集venn文件", 0],
            [r".*venn.*\.pdf", "", "差异基因集venn图", 0],
            [r'.*upset.*\.pdf', 'txt', '差异基因集upset图', 0]
        ])
        super(DiffGenesetVennWorkflow, self).end()

    def set_db(self):
        all_exp = self.api.api("medical_transcriptome.diff_geneset")
        all_exp.add_venn_tt(
            self.option('main_id'),
        )

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        chart.chart_geneset_venn(self.option("geneset_id"))
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))