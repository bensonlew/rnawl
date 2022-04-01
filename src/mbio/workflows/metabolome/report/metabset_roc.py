# -*- coding: utf-8 -*-

#zouguanqing
""""""

import datetime
from biocluster.workflow import Workflow
import re
import os
import json
import shutil
import pandas as pd

class MetabsetRocWorkflow(Workflow):
    """
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetRocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "exp_table", "type": "infile", "format": "sequence.profile_table"},  #表达量表
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组表，只接受两个分组。接口用to_file 生成
            {"name": "group_detail", "type": "string"},
            {"name": "roc_id", "type": "string"},  # 接口传主表id
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            #{"name": "compound", "type": "string"},  #选择代谢物名称
            {"name": "update_info", "type": "string"},
            {"name": "confidence_interval", "type": "float", "default": 0.95},
            {'name': 'save_pdf', 'type': 'int', 'default': 0},  # 是否保存pdf图片导结果目录，v3.5升级
            {"name": "name", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.option("save_pdf", 0)
        self.roc = self.add_module("metabolome.roc_batch")
        self.output_dir = self.roc.output_dir

    def run_roc(self):
        options = {
            'exp_table': self.option('exp_table'),
            'group_table': self.option('group_table'),
            'metab_set_table' : self.option('metab_set_table'),
            'confidence_interval': str(self.option('confidence_interval')),
            'metab_desc' : self.option('metab_desc')
        }
        self.roc.set_options(options)
        self.roc.on('end', self.set_db)
        self.output_dir = self.roc.output_dir
        self.roc.run()

    def end(self):
        if self.option("save_pdf"):
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.roc.output_dir))
        result_dir = self.add_upload_dir(self.roc.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ROC分析结果目录", 0],
            ["./roc_curve.xls", "xls", "ROC曲线结果表", 0],
            ["./roc_auc.xls", "xls", "AUC值计算结果表", 0],
            ["./roc_interval.xls", "xls", "ROC曲线置信区间", 0],
            ["./best_loc.xls", "xls", "ROC曲线临界值", 0],
            ["./roc_curve_smooth.xls", "xls", "ROC曲线平滑处理结果表", 0]
        ])
        super(MetabsetRocWorkflow, self).end()

    def set_db(self):
        api_roc = self.api.api('metabolome.metabset_roc')
        datacurve = self.output_dir + '/roc_curve.xls'
        api_roc.add_roc_curve(roc_id=self.option("roc_id"), type="", file=datacurve)

        # if os.path.exists(self.output_dir + '/roc_curve_smooth.xls'):
        #     datacurve_s = self.output_dir + '/roc_curve_smooth.xls'
        #     api_roc.add_roc_curve(roc_id=self.option("roc_id"), type="smooth", file=datacurve_s)

        if os.path.exists(self.output_dir + '/roc_interval.xls'):
            datainterval = self.output_dir + '/roc_interval.xls'
            api_roc.add_roc_interval(roc_id=self.option("roc_id"), file=datainterval)

        dataauc = self.output_dir + '/roc_auc.xls'
        api_roc.add_roc_auc(roc_id=self.option("roc_id"), file=dataauc, type="")

        # if os.path.exists(self.output_dir + '/roc_auc_smooth.xls'):
        #     dataauc_s = self.output_dir + '/roc_auc_smooth.xls'
        #     api_roc.add_roc_auc(roc_id=self.option("roc_id"), file=dataauc_s, type="smooth")
        databestloc = self.output_dir + '/best_loc.xls'
        api_roc.add_roc_best_loc(roc_id=self.option("roc_id"), file=databestloc)
        #v3 add insert metabset_roc_table
        api_roc.add_roc_table(self.option("roc_id"), metab_desc=self.option('metab_desc').path)

        report_files = [datacurve, dataauc, databestloc]
        for f in report_files:
            if not os.path.isfile(f):
                self.logger.error("找不到报告文件:{}".format(f))
                self.set_error("找不到报告文件")
        # self.option("save_pdf", 0)
        if self.option("save_pdf"):
            self.figsave = self.add_tool("metabolome.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("roc_id"),
                "table_name": self.option("name"),
                "project": "metabolome",
                "submit_loc": "metabsetroc",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        self.run_roc()
        super(MetabsetRocWorkflow, self).run()


# if __name__ == '__main__':
#     from biocluster.sheet import Sheet
#
#     data = {
#         "name" : "roc",
#         "id" : "tsg_36964",
#         'type': 'workflow',
#         "options" : {
#             "exp_table" : ,
#             "metab_desc" : ,
#             "group_table" : ,
#             "metab_set_table" : ,
#             "roc_id" : ""
#         }
#
#     }
#
#     data = Sheet(data=data)
#     wf = MetabsetRocWorkflow(data)
#     wf.run()