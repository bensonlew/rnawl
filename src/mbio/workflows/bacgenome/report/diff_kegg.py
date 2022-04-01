# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# __modify__ = '2019/4'

from biocluster.workflow import Workflow
import os
import shutil

class DiffKeggWorkflow(Workflow):
    """
    kegg蛋白差异分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "samples", "type": "string"},
            {"name": "pathway_id", "type": "string"},
            {"name": "main_table_id","type":"string"},
            {"name": "k_list","type":"string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.diff_kegg = self.add_tool("bacgenome.diff_kegg")

    def run_diff_kegg(self):
        self.logger.info("start run kegg!")
        options = {
            "pathway_id": self.option('pathway_id'),
            "k_list": self.option('k_list')
        }
        self.diff_kegg.set_options(options)
        self.diff_kegg.on('end', self.set_db)
        self.diff_kegg.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run diff kegg workflow")
        self.run_diff_kegg()
        super(DiffKeggWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        if os.path.exists(self.output_dir + "/pathway_img"):
            shutil.rmtree(self.output_dir + "/pathway_img")
        shutil.copytree(self.diff_kegg.output_dir+'/pathway_img', self.output_dir + "/pathway_img")
        self.logger.info("开始进行导表")
        kegg_api = self.api.api('bacgenome.diff_kegg')
        kegg_api.add_diff_kegg_detail(self.option("main_table_id"),self.option("samples"), self.output_dir + "/pathway_img")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢通路差异分析结果目录",0],
            ["pathway_img", "", "通路图片",0],
        ])
        super(DiffKeggWorkflow, self).end()




