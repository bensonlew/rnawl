# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import os
import pandas as pd
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.prok_rna.chart import Chart


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
            dict(name="pca_main_id", type='string'),
            dict(name="exp_level", type='string'),
            dict(name="group_table", type="string"),
            dict(name="analysis_type", type="string", default="pca"),
            # dict(name="type", type='string'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.exp_pca")
        self.ellipse = None
        group_dict = json.loads(self.option("group_dict"))
        i = 0
        for sample_infos in group_dict.values():
            if len(sample_infos) < 3:
                i += 1
            else:
                continue
        if i == 0:
            self.ellipse = self.add_tool("graph.ellipse")

    def run(self):
        # self.tool.on("end", self.set_db)
        # self.run_tool()
        if not self.ellipse is None:
            self.tool.on("end", self.run_ellipse)
            self.ellipse.on("end", self.set_db)
            self.run_tool()
        else:
            self.tool.on("end", self.set_db)
            self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def run_ellipse(self):
        options = {}
        if self.option("group_table"):
            options['group_table'] = self.option("group_table")
        pc_map = {'pca':"/PCA.xls",'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda':'/Dbrda/db_rda_sites.xls','nmds':'/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  ##'rda_cca'
                  }
        options['analysis'] = self.option('analysis_type')
        options['pc_table'] = self.tool.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(options)
        self.ellipse.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("prok_rna.all_exp")
        # add result info
        all_exp.add_exp_pca2(self.tool.work_dir, main_id=self.option('pca_main_id'), exp_level=self.option('exp_level'))
        if not self.ellipse is None:
           all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('pca_main_id')))
        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        group_dict = json.loads(self.option("group_dict"))
        exp_pca_file = os.path.join(self.tool.output_dir, 'PCA.xls')
        exp_pca_var_file = os.path.join(self.tool.output_dir, 'Explained_variance_ratio.xls')
        if self.ellipse:
            exp_pca_ellipse = os.path.join(self.ellipse.output_dir, 'ellipse_out.xls')
        else:
            exp_pca_ellipse = None
        chart.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict, exp_pca_ellipse=exp_pca_ellipse, pcs=["PC1", "PC2"])
        chart.to_pdf()

        # move pdf
        pca_pdf = os.path.join(self.work_dir, 'all.exp_relation_pca.scatter.pdf')
        new_path = os.path.join(self.tool.output_dir, 'sample_pca.pdf')
        self.move_pdf(pca_pdf, new_path)
        if os.path.exists(os.path.join(self.work_dir, 'all.exp_relation_pca_ell.scatter.pdf')):
            new_path = os.path.join(self.tool.output_dir, 'sample_pca_ell.pdf')
            self.move_pdf(os.path.join(self.work_dir, 'all.exp_relation_pca_ell.scatter.pdf'), new_path)

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.tool.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PCA分析结果目录"],
            ["PCA.xls", "", "样本间pca统计表"],
            ["Explained_variance_ratio.xls ", "", "样本间pca解释率统计表"],
            ["sample_pca.pdf", "pdf", "样本间PCA图"],
            ["sample_pca_ell.pdf", "pdf", "样本间PCA图（置信圈）"],
        ])
        super(ExpPcaWorkflow, self).end()

    def run_tool(self):
        options = dict(
            exp=self.option('exp_matrix'),
        )
        self.tool.set_options(options)
        self.tool.run()
