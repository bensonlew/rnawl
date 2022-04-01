# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
from mbio.packages.prok_rna.chart import Chart
import glob


class ExpVennWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpVennWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="venn_main_id", type='string'),
            dict(name="threshold", type='string'),
            dict(name="exp_level", type='string'),
            dict(name="group", type="string"),
            # dict(name="type", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.exp_venn")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(ExpVennWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        all_exp = self.api.api("prok_rna.all_exp")
        # add result info
        graph_table = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        all_exp.add_exp_venn(graph_table, main_id=self.option('venn_main_id'), exp_level=self.option('exp_level'))
        self.end()

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        exp_venn = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        chart.chart_exp_venn(exp_venn)
        chart.to_pdf()

        # move pdf
        venn_pdf = os.path.join(self.work_dir, 'all.exp.venn.pdf')
        new_path = os.path.join(self.output_dir, 'sample_venn.pdf')
        self.move_pdf(venn_pdf, new_path)

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "样本间Venn分析结果目录"],
            ['sample_venn.pdf', "", "样本间Venn图"],
        ])
        super(ExpVennWorkflow, self).end()

    def run_tool(self):
        options = dict(
            express_matrix=self.option('exp_matrix'),
            group_table=self.option('group'),
            threshold=self.option('threshold'),
        )
        self.tool.set_options(options)
        self.tool.run()
