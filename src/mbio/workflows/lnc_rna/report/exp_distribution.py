# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time


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
            dict(name="type", type="string"),
            dict(name="rna_type",type="string"),
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

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # time.sleep(1)
        all_exp = self.api.api("lnc_rna.all_exp")
        all_exp.add_distribution(
            self.option('exp_matrix'),
            main_id=self.option('graph_main_id'),
            group_dict=json.loads(self.option("group_dict")),
        )

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "表达量分布分析结果目录"],
        # ])
        super(ExpDistributionWorkflow, self).end()

