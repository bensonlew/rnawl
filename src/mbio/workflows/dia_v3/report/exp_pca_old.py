# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import pandas as pd
import os

class ExpPcaWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="group", type="string"),
            dict(name="pca_main_id", type='string'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("labelfree.pca")


    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("dia.pca")
        # add result info
        all_exp.add_pca(self.tool.output_dir, main_id=self.option('pca_main_id'),)
        self.end()

    def end(self):
        os.remove(self.tool.output_dir + "/pca_rotation.xls")
        os.rename(self.tool.output_dir + "/pca_rotation_all.xls", self.tool.output_dir + "/pca_rotation.xls")
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
        [".", "", "样本PCA分析结果目录"],
        ["./pca_importance.xls", "xls", "PCA主成分解释表"],
        ["./pca_rotation.xls", "xls", "蛋白相关主成分贡献度"],
        ["./pca_sites.xls", "xls", "PCA分析样本坐标表"],
        ])
        super(ExpPcaWorkflow, self).end()

    def run_tool(self):
        options = dict(
            otutable=self.option('exp_matrix'),
        )
        self.tool.set_options(options)
        self.tool.run()
