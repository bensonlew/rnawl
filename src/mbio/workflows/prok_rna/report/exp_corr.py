# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import os
import pandas as pd
import unittest
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.prok_rna.chart import Chart
from collections import OrderedDict


class ExpCorrWorkflow(Workflow):
    """
    表达量相关性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name='scm', type='string', default='complete'),
            dict(name='scd', type='string', default='correlation'),
            dict(name='corr_method', type='string', default='pearson'),
            dict(name="corr_main_id", type='string'),
            dict(name="exp_level", type='string'),
            dict(name='log_base', type='string', default='no'),
            dict(name='Draw_in_groups', type='string', default='no'),
            # dict(name="type", type='string'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.exp_corr")
        if self.option("Draw_in_groups") == "yes":
            group_dict = json.loads(self.option("group_dict"))
            ori_exp=pd.read_table(self.option('exp_matrix'), index_col="seq_id",header=0)
            for i in group_dict.keys():
                t = i + "=" + "(" + group_dict[i][0]
                if len(group_dict[i]) > 1:
                    for j in group_dict[i][1:]:
                        t = t + "+" + j
                t = t + ")" + "/" + str(len(group_dict[i]))
                ori_exp.eval(t, inplace=True)
            final_exp = ori_exp.loc[:, group_dict.keys()]
            final_exp.to_csv(self.work_dir+"/work_exp_matrix", sep='\t')
            self.final_exp = self.work_dir+"/work_exp_matrix"
        else:
            self.final_exp = self.option('exp_matrix')

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(ExpCorrWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("prok_rna.all_exp")
        # add result info
        all_exp.add_exp_corr2(self.tool.work_dir, main_id=self.option('corr_main_id'), exp_level=self.option('exp_level'))
        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        exp_corr_file = os.path.join(self.tool.work_dir, 'sample_correlation.xls')
        exp_corr_tree_file = os.path.join(self.tool.work_dir, 'sample.cluster_tree.txt')
        if self.option("Draw_in_groups") == "yes":
            group_dict = json.loads(self.option("group_dict"))
            group_dict_chart = OrderedDict([(i, [i]) for i in group_dict.keys()])
            group_dict_chart = json.dumps(group_dict_chart)
        else:
            group_dict_chart = self.option('group_dict')
        chart.chart_exp_corr(exp_corr_file, exp_corr_tree_file, group_dict_chart)
        chart.to_pdf()

        # move pdf
        self.move_pdf(os.path.join(self.work_dir, 'exp.heatmap.heat_corr.pdf'),
                      os.path.join(self.tool.output_dir, 'sample_correlation.pdf'))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性分析结果目录"],
            ["sample_correlation.xls", "", "样本间相关性系数统计表"],
            ["expression_matrix.xls", "", "样本间相关性表达矩阵"],
            ["sample_correlation.pdf", "pdf", "样本间相关性热图"],
        ])
        super(ExpCorrWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.final_exp,
            scm=self.option('scm'),
            scd=self.option('scd'),
            corr_method=self.option('corr_method'),
        )
        if self.option('log_base') == 'no':
            pass
        else:
            options.update({
                'log_base': int(self.option('log_base')) 
            })
        self.tool.set_options(options)
        self.tool.run()
