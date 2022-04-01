# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
from mbio.packages.prok_rna.chart import Chart
import os
import glob
from collections import OrderedDict


class ExpDistributionWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpDistributionWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="graph_main_id", type="string"),
            dict(name="group_dict", type="string"),
            dict(name="exp_matrix", type="string"),
            dict(name="exp_level", type="string"),
            # dict(name="type", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        # 没有tool时，意味着不需要self.run, 进一步意味着需要发起监听。
        self.start_listener()
        self.fire("start")
        self.set_db()
        self.end()
        # super(ExpDistributionWorkflow, self).run() 不需要

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        group_dict = json.loads(self.option("group_dict"), object_pairs_hook=OrderedDict)
        samples = list()
        for each in group_dict.keys():
            samples.extend(group_dict[each])
        chart.chart_exp_dis(self.option('exp_matrix'), group_dict, samples)
        chart.to_pdf()

        # move pdfs
        distribution_pdfs = glob.glob(os.path.join(self.work_dir, "*exp_distribution*.pdf"))
        for each in distribution_pdfs:
            dis_type = os.path.basename(each).split('.')[0]
            if each.endswith('.density.pdf'):
                new_file = dis_type + '_exp_density.pdf'
            elif each.endswith('.box.pdf'):
                new_file = dis_type + '_exp_box.pdf'
            elif each.endswith('.violin.pdf'):
                new_file = dis_type + '_exp_violin.pdf'
            new_path = os.path.join(self.output_dir, new_file)
            self.move_pdf(each, new_path)

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # time.sleep(1)
        all_exp = self.api.api("prok_rna.all_exp")
        all_exp.add_distribution(
            self.option('exp_matrix'),
            main_id=self.option('graph_main_id'),
            exp_level=self.option('exp_level'),
            group_dict=json.loads(self.option("group_dict")),
        )

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达量分布分析结果目录"],
        ])
        result_dir.add_regexp_rules([
            [r".*_exp_density\.pdf", "pdf", "表达量分布密度图"],
            [r".*_exp_box\.pdf", "pdf", "表达量分布盒形图"],
            [r".*_exp_violin\.pdf", "pdf", "表达量分布小提琴图"],
        ])
        super(ExpDistributionWorkflow, self).end()

