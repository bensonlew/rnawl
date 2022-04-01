# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest
import re
from mbio.packages.whole_transcriptome.utils import check_map_dict
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from biocluster.workflow import Workflow
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.whole_transcriptome.chart.chart import Chart
import glob
from shutil import copyfile

class ExpPcaWorkflow(Workflow):
    '''
    last_modify: 2019.11.14
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        LIBRARY = ('long', 'small', 'circle')
        LEVEL = ('G', 'T')
        KIND = ('all', 'ref', 'new')
        options = [
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'library', 'type': 'string', 'default': LIBRARY[0]},
            {'name': 'level', 'type': 'string', 'default': LEVEL[1]},
            {'name': 'kind', 'type': 'string', 'default': KIND[0]},
            {'name': 'group_id', 'type': 'string', 'default': None},
            {'name': 'group_dict', 'type': 'string', 'default': None},
            {'name': 'take_mean', 'type': 'bool', 'default': False},
            {'name': 'is_rmbe', 'type': 'string', 'default': 'false'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.inter_dirs = []

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(ExpPcaWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
        self.run_exp_build()
        super(ExpPcaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="exp_pca", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_exp_build(self):
        self.exp_build = self.add_tool('whole_transcriptome_v1_1.formation.exp_build')
        self.exp_build.set_options({
            'task_id': self.option('task_id'),
            'library': self.option('library'),
            'level': self.option('level'),
            'kind': self.option('kind'),
            'group_dict': self.option('group_dict'),
            'is_rmbe': self.option('is_rmbe')
        })
        self.exp_build.on('end', self.run_exp_pca)
        self.exp_build.run()

    def run_exp_pca(self):
        self.exp_pca = self.add_tool('whole_transcriptome.formation.exp_pca')
        self.exp_pca.set_options({
            'exp_matrix': self.exp_build.option('exp_matrix'),
            'group_table': self.exp_build.option('group_table').path,
            'take_mean': self.option('take_mean')
        })
        self.exp_pca.on('end', self.set_output)
        self.exp_pca.run()

    def set_output(self):
        for fname in os.listdir(self.exp_pca.output_dir):
            source = os.path.join(self.exp_pca.output_dir, fname)
            link_name = os.path.join(self.output_dir, fname)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.set_db()

    def set_db(self):
        api = self.api.api('whole_transcriptome.expression')
        map_dict = {
            'evr': os.path.join(self.output_dir, 'explained_variance_ratio.txt'),
            'pca': os.path.join(self.output_dir, 'pca.txt'),
            'ellipse': os.path.join(self.output_dir, 'ellipse.txt')
        }
        if not os.path.isfile(map_dict['ellipse']):
            map_dict.pop('ellipse')
        task_id = self.option('task_id')
        library = self.option('library')
        level = self.option('level')
        kind = self.option('kind')
        group_id = self.option('group_id')
        group_dict = json.loads(self.option('group_dict'))
        take_mean = 'yes' if self.option('take_mean') else 'no'
        pca_id = self.option('main_id')
        way = api.db['exp'].find_one({'task_id': task_id, 'level': level})['way']
        project_sn = self.sheet.project_sn
        api.add_exp_pca(map_dict, library, level, kind, group_id, group_dict, way, take_mean, task_id, project_sn,
                        pca_id)
        self.set_upload()

    def Chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        group_dict = self.option('group_dict')
        group_dict = json.loads(self.option('group_dict'))
        lib_type = self.option('library')
        exp_pca_file = os.path.join(self.output_dir, 'pca.txt')
        exp_pca_var_file = os.path.join(self.output_dir, 'explained_variance_ratio.txt')
        exp_pca_ellipse = os.path.join(self.output_dir, 'ellipse.txt')
        if not os.path.isfile(exp_pca_ellipse):
            exp_pca_ellipse = None
        chart.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict,
                           exp_pca_ellipse=exp_pca_ellipse, pcs=["PC1", "PC2"], lib_type=lib_type)
        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        upload_dir = os.path.join(self.work_dir, 'upload')
        for p in pdf_file:
            copyfile(p, os.path.join(upload_dir, os.path.basename(p)))

    def set_upload(self):

        if self.option('library') == 'long':
            from mbio.packages.whole_transcriptome.catalogue import mrna as rna
        elif self.option('library') == 'small':
            from mbio.packages.whole_transcriptome.catalogue import mirna as rna
        elif self.option('library') == 'circle':
            from mbio.packages.whole_transcriptome.catalogue import circrna as rna
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/01 Express/02 Exp_PCA')
        upload_dir = os.path.join(self.work_dir, 'upload')
        if not os.path.isdir(upload_dir):
            os.mkdir(upload_dir)
        self.Chart()
        rna.set_exp_pca(self.option('task_id'), upload_dir)
        if os.path.exists(os.path.join(upload_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(upload_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(upload_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(upload_dir)
        self.inter_dirs = [
            ["01 Express", "", "表达量结果目录",0],
            ["01 Express/02 Exp_PCA", "", "样本间PCA分析", 0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '样本间PCA分析文件', 0],
            ['explained_variance_ratio.xls', 'xls', '主成分解释表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['*pca.*pdf', 'pdf', '样本间PCA图', 0],
            ['*pca_ell*.pdf', 'pdf', '样本间PCA带置信圈图', 0],
        ])
        self.end()

    def end(self):
        super(ExpPcaWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.whole_transcriptome.report.exp_pca import ExpPcaWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'exp_pca_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.report.exp_pca',
            'options': {
                'task_id': 'tsg_36088',
                'library': 'long',
                'level': 'T',
                'kind': 'new',
                'group_id': '5dcb545a17b2bf08b8f25ced',
                'group_dict': json.dumps({'Acute': ['A1', 'A2'], 'Stable': ['S1', 'S3'], 'control': ['C2', 'C3']}),
                'take_mean': 'no',
            }
        }
        wsheet = Sheet(data=data)
        wf = ExpPcaWorkflow(wsheet)
        wf.sheet.project_sn = '188_5dba6f542345b'
        wf.sheet.task_id = 'tsg_36088'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
