# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.workflow import Workflow

from mbio.packages.whole_transcriptome.utils import check_map_dict
from mbio.packages.whole_transcriptome.chart.chart import Chart
from biocluster.config import Config
from bson import ObjectId
import glob
from shutil import copyfile

class ExpGraphWorkflow(Workflow):
    '''
    last_modify: 2019.11.13
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpGraphWorkflow, self).__init__(wsheet_object)
        LEVEL = ('G', 'T')
        CATEGORY = ('mRNA', 'lncRNA', 'miRNA', 'circRNA')
        KIND = ('all', 'ref', 'new')
        options = [
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': LEVEL[1]},
            {'name': 'category', 'type': 'string', 'default': CATEGORY[0]},
            {'name': 'kind', 'type': 'string', 'default': KIND[0]},
            {'name': 'group_id', 'type': 'string', 'default': None},
            {'name': 'group_dict', 'type': 'string', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def run(self):
        self.run_exp_build()
        super(ExpGraphWorkflow, self).run()

    def run_exp_build(self):
        self.exp_build = self.add_tool('whole_transcriptome.formation.exp_build')
        self.exp_build.set_options({
            'task_id': self.option('task_id'),
            'level': self.option('level'),
            'category': self.option('category'),
            'kind': self.option('kind'),
            'group_dict': self.option('group_dict')
        })
        self.exp_build.on('end', self.run_exp_graph)
        self.exp_build.run()

    def run_exp_graph(self):
        self.exp_graph = self.add_tool('whole_transcriptome.formation.exp_graph')
        self.exp_graph.set_options({
            'exp_matrix': self.exp_build.option('exp_matrix'),
            'group_table': self.exp_build.option('group_table')
        })
        self.exp_graph.on('end', self.set_output)
        self.exp_graph.run()

    def set_output(self):
        for fname in os.listdir(self.exp_graph.output_dir):
            source = os.path.join(self.exp_graph.output_dir, fname)
            link_name = os.path.join(self.output_dir, fname)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.set_db()

    def set_db(self):
        api = self.api.api('whole_transcriptome.expression')
        map_dict = {
            'sample_box': os.path.join(self.output_dir, 'sample_box_data.pk'),
            'group_box': os.path.join(self.output_dir, 'group_box_data.pk'),
            'sample_density': os.path.join(self.output_dir, 'sample_density_data.pk'),
            'group_density': os.path.join(self.output_dir, 'group_density_data.pk'),
            'sample_volin': os.path.join(self.output_dir, 'sample_volin_data.pk'),
            'group_volin': os.path.join(self.output_dir, 'group_volin_data.pk'),
        }
        # check_map_dict(map_dict)
        task_id = self.option('task_id')
        level = self.option('level')
        category = self.option('category')
        kind = self.option('kind')
        group_id = self.option('group_id')
        group_dict = json.loads(self.option('group_dict'))
        graph_id = self.option('main_id')
        exp_id = api.db['exp'].find_one({'task_id': task_id, 'level': level})['main_id']
        project_sn = self.sheet.project_sn
        api.add_exp_graph(map_dict, category, exp_id, level, kind, group_id, group_dict, task_id, project_sn, graph_id)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        group_dict = json.loads(self.option("group_dict"))
        exp = self.option("exp_matrix")

        # 获取exp_type
        exp_type = "TPM"
        print
        "exp_type", exp_type
        # 绘制分组均值分布图
        chart.chart_exp_dis_one(exp, group_dict, exp_type)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            if os.path.exists(
                    self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "group_distribution")):
                os.remove(self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "group_distribution"))
            copyfile(p, self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "group_distribution"))
        # 绘制样本分布图

        chart.chart_exp_dis_one(exp, group_dict=None, exp_type=exp_type)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            if os.path.exists(
                    self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution")):
                os.remove(
                    self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))
            copyfile(p, self.output_dir + "/" + os.path.basename(p).replace("exp_distribution", "sample_distribution"))

    def end(self):
        self.chart()
        super(ExpGraphWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.whole_transcriptome.report.exp_graph import ExpGraphWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'exp_graph_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.report.exp_graph',
            'options': {
                'task_id': 'tsg_36088',
                'level': 'T',
                'category': 'mRNA',
                'kind': 'all',
                'group_id': '5dcb545a17b2bf08b8f25ced',
                'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}),
            }
        }
        wsheet = Sheet(data=data)
        wf = ExpGraphWorkflow(wsheet)
        wf.sheet.project_sn = '188_5dba6f542345b'
        wf.sheet.task_id = 'tsg_36088'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
