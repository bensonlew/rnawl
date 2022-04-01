# -*- coding: utf-8 -*-
import shutil
import os
import gevent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.common import clean_mongo

class DeleteRelationWorkflow(Workflow):
    """
    1.用于页面删除交互分析的记录同时删除mongo数据；
    2.用于页面删除task或者批量删除，同时删除mongo数据；
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeleteRelationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},#删除任务的id
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check(self):
        if not self.option('task_id'):
            raise OptionError("Must set task id of what you want to delate")
        if not self.option('project_type'):
            raise OptionError("Must set project type of what you want to delate")

    def delete_relation(self):
        """删除所有关联分析记录"""
        relation_tables = ["relation_corr", "relation_ipath", "relation_keggp", "relation_keggpview", "relation_kegg_heatmap",
                           "relation_corr_network", "relation_o2pls", "relation_procrustes", "relation_keggp_enrich"]
        clean_mongo({"task_id": self.option("task_id")}, relation_tables)
        gevent.spawn_later(5, self.end)

    def run(self):
        self.delete_relation()
        super(DeleteRelationWorkflow, self).run()

    def end(self):
        super(DeleteRelationWorkflow, self).end()