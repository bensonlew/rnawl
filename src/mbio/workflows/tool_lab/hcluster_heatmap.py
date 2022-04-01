# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""cluster_heatmap"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import pandas as pd
from biocluster.file import download,exists
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.common_function import meta_get_file
import shutil

class HclusterHeatmapWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HclusterHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":'otutable', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {"name": "sep", "type": 'string', "default": "tab"},
            {"name": "group", "type": 'infile', "format": "ref_rna_v2.common"},
            {"name": 'log_change', "type": 'string'},
            # {"name": '0_mean', 'type': 'bool'},
            # {"name": 'zscore', 'type': 'bool'},
            {'name': 'scale', 'type': 'string', },
            {"name": 'cluster', 'type': 'bool'},
            {'name': 'gcd', 'type': 'string'},   #行距离算法
            {"name": 'gcm', 'type': 'string'},   #行聚类方式
            {"name": 'scd', 'type': 'string'},   #列距离算法
            {"name": 'scm', 'type': 'string'},   #列聚类方式
            {"name": 'title', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'project_data', 'type': 'string', 'default': ''},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'project_task_id', 'type': 'string'},

        ]
        self.add_option(options)
        self.revise_infiles()
        # self.tool = self.add_tool("tool_lab.table_standard")
        self.table_select = self.add_tool("tool_lab.hcluster_heatmap.table_select")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option('project_data'):
            self.data_table,self.group_table = meta_get_file(self.option('project_data'),self.work_dir)
        else:
            self.data_table = self.option('otutable').path
            self.group_table = self.option('group').path
        self.file_check()
        self.run_table_select()
        super(HclusterHeatmapWorkflow, self).run()

    # added by zhangyitong on 20210909
    def file_check(self):
        with open(self.data_table, 'r') as f:
            cols = f.readline().strip().split('\t')
        if len(cols) > len(set(cols)):
            self.set_error('输入数据表含有重复样本名（列名），请检查。')

        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
        in_df = pd.read_table(self.data_table, sep=sep)
        num_row, num_col = in_df.shape
        if self.option('cluster') and num_col > 3000:
            self.set_error('输入数据样本（列）超过3000个，暂不支持层级聚类。')
        if self.option('cluster') and num_row > 3000:
            self.set_error('输入数据特征（行）超过3000个，暂不支持层级聚类。')

    def run_table_select(self):
        if self.option('sep').lower() == 'blank':
            sep = 'space'
        else:
            sep = self.option('sep').lower()
        options = {
            "origin_table": self.data_table,
            "sep": sep,
        }
        self.table_select.set_options(options)
        self.table_select.on('end', self.run_hcluster_heatmap)
        self.table_select.run()

    def run_hcluster_heatmap(self):
        if self.option('sep').lower() == 'blank':
            sep = 'space'
        else:
            sep = self.option('sep').lower()
        self.hcluster_heatmap = self.add_module('tool_lab.hcluster_heatmap')
        if self.option('cluster') == False:
            self.hcluster_heatmap.set_options({
                "otutable": self.table_select.option('select_table'),
                "sep": sep,
                "log_change": self.option('log_change'),
                'scale': self.option('scale'),
                'cluster': self.option('cluster'),
            })
        else:
            self.hcluster_heatmap.set_options({
                "otutable": self.table_select.option('select_table'),
                "sep": sep,
                "log_change": self.option('log_change'),
                'scale': self.option('scale'),
                'cluster': self.option('cluster'),
                'scd': self.option('scd'),
                'scm': self.option('scm'),
                'gcd': self.option('gcd'),
                'gcm': self.option('gcm')
            })
        self.hcluster_heatmap.on('end', self.set_db)
        self.hcluster_heatmap.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        hcluster = self.api.api("tool_lab.hcluster_heatmap")
        # add result info
        sample_tre = os.path.join(self.hcluster_heatmap.output_dir, 'sample_hcluster.tre')
        feature_tre = os.path.join(self.hcluster_heatmap.output_dir, 'feature_hcluster.tre')
        heatmap = os.path.join(self.hcluster_heatmap.output_dir, 'scale_table.txt')
        hcluster.add_hcluster_heatmap(sample_tre, feature_tre, heatmap, self.group_table, main_id=self.option('main_id'), title=self.option('title'))
        self.end()

    def end(self):
        if self.option("project_data"):
            if os.path.exists(self.hcluster_heatmap.output_dir + "/input_data"):
                shutil.rmtree(self.hcluster_heatmap.output_dir + "/input_data")
            os.mkdir(self.hcluster_heatmap.output_dir + "/input_data")
            os.link(self.data_table, self.hcluster_heatmap.output_dir + "/input_data/input_table.txt")
            os.link(self.group_table, self.hcluster_heatmap.output_dir + "/input_data/input_group.txt")
        result_dir = self.add_upload_dir(self.hcluster_heatmap.output_dir)
        super(HclusterHeatmapWorkflow, self).end()
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.hcluster_heatmap import HclusterHeatmapWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'cluster_heatmap_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.hcluster_heatmap',
            'options': {
                "otutable" : "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/exp_matrix_new",
                'group':"/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/group_1.txt",
                "sep": 'tab',
                "log_change": None,
                'scale': '0_mean',
                'cluster': True,
                'scd': 'euclidean',
                'scm': 'average',
                'gcd': 'euclidean',
                'gcm': 'average',
                'title': 'test_cluster_heatmap'
            }
        }
        wsheet = Sheet(data=data)
        wf =HclusterHeatmapWorkflow(wsheet)
        wf.sheet.id = 'cluster_heatmap'
        wf.sheet.project_sn = 'cluster_heatmap'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
