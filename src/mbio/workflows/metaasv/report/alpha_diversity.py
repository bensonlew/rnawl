# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.workflow import Workflow
import os
from mbio.api.to_file.meta import *
from mbio.packages.metaasv.filter_newick import get_level_newicktree
from mbio.packages.metaasv.common_function import link_file
import re

class AlphaDiversityWorkflow(Workflow):
    """
    metaasv Alpha多样性指数
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        # self.rpc = False   AlphaDiversityWorkflow
        super(AlphaDiversityWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的OTU id
            {"name": "asv_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "indices", "type": "string", "default": "sobs"},
            {"name": "level", "type": "int","default": 9},
            {"name": "main_id", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},#输入Group表
            {"name": "group_detail", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.estimators = self.add_tool('meta.alpha_diversity.estimators')
        self.indices = self.option("indices")

    def run(self):
        """
        运行
        :return:
        """
        self.get_tree()
        super(AlphaDiversityWorkflow, self).run()

    def run_estimators(self):
        indices = self.option('indices').split(',')
        if 'pd' in indices:
            ####如果需要得到PD指数，则需要找到OTU对应的进化树
            opts = {
                'otu_table': self.tree.option("temp_otu_file").prop["path"],
                'newicktree': self.tree.option("temp_tree_file").prop["path"],
                'indices': self.option('indices')
            }
        else:
            opts = {
                'otu_table': self.option('otu_file'),
                'indices': self.option('indices')
            }
        self.logger.info(self.option('indices'))
        self.estimators.set_options(opts)
        self.estimators.on('end', self.set_db)
        self.estimators.run()
        self.output_dir = self.estimators.output_dir

    def get_tree(self):
        indices = self.option('indices').split(',')
        if 'pd' in indices:
            self.tree = self.add_tool('metaasv.get_level_newicktree')
            opts = {
                'otu_file': self.option('otu_file').prop["path"],
                'asv_id': self.option('asv_id'),
                'level': self.option('level')
            }
            self.tree.set_options(opts)
            self.tree.on('end', self.run_estimators)
            self.tree.run()
        else:
            self.run_estimators()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        # link_dir(self.estimators.output_dir, self.output_dir)
        api_estimators = self.api.api("metaasv.estimator")
        est_path = self.output_dir + "/estimators.xls"
        link_file(os.path.join(self.estimators.work_dir,"estimators.xls"), est_path)
        if not os.path.isfile(est_path):
            self.logger.error("找不到报告文件:{}".format(est_path))
            self.set_error("找不到报告文件")
        main_id = api_estimators.add_est_table(est_path, level=self.option('level'),otu_id=self.option('asv_id'),
                                              est_id=self.option("main_id"),indices=self.option("indices"))
        api_estimators.add_est_bar(est_path, self.option("main_id"),indices=self.option("indices"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "多样性指数结果目录", 0, ""],
            ["./estimators.xls", "xls", "alpha多样性指数表", 0, ""]
        ])
        # print self.get_upload_files()
        super(AlphaDiversityWorkflow, self).end()
