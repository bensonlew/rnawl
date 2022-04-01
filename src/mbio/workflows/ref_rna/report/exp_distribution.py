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
            # dict(name="group_dict", type="string"),
            dict(name="exp_matrix", type="string"),
            dict(name="express_level", type="string"),
            dict(name="type", type="string"),
            dict(name="group_id", type="string", default=""),
            dict(name="group_detail", type="string", default=""),
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

    def get_group_dict(self):
        group_file = self.option("group_id")
        if not group_file:
            return None
        group_dict = dict()
        group_file = group_file.split(',')[0]
        with open(group_file) as f:
            _ = f.readline()
            for line in f:
                sample, group = line.strip().split()
                group_dict.setdefault(group, list())
                group_dict[group].append(sample)
        return group_dict

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        # time.sleep(1)
        all_exp = self.api.api("ref_rna.distribution")
        group_dict = self.get_group_dict()
        exp_matrix = self.option("exp_matrix").split(',')[0]
        all_exp.add_distribution(
            exp_matrix,
            main_id=self.option('graph_main_id'),
            group_dict=group_dict,
        )

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "表达量分布分析结果目录"],
        # ])
        super(ExpDistributionWorkflow, self).end()

