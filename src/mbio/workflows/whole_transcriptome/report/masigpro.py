# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

from biocluster.workflow import Workflow
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class MasigproWorkflow(Workflow):
    '''
    last_modify: 2019.12.18
    '''

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MasigproWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.matrix'},
            {'name': 'geneset', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'design', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'cluster', 'type': 'int', 'default': 9},
            {'name': 'method', 'type': 'string', 'default': 'hclust'},  # ['hclust', 'kmeans', 'Mclust']
            {'name': 'level', 'type': 'string', 'default': 'T'},  # ['T', 'G']
            {'name': 'way', 'type': 'string', 'default': 'tpm'},  # ['tpm', 'fpkm']
            {'name': 'group_dict', 'type': 'string'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/05 Time_series_Analysis')
        self.inter_dirs = []

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
        super(MasigproWorkflow, self).send_log(data)

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.get_run_log()
        self.run_masigpro()
        super(MasigproWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="masigpro",
                                main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_masigpro(self):
        self.step.add_steps('masigpro')
        self.masigpro = self.add_tool('whole_transcriptome.timeseries.masigpro')
        options = {
            'matrix': self.option('matrix'),
            'geneset': self.option('geneset'),
            'design': self.option('design'),
            'cluster': self.option('cluster'),
            'method': self.option('method'),
            'level': self.option('level'),
            'way': self.option('way'),
        }
        self.masigpro.set_options(options)
        self.masigpro.on('start', self.set_step, {'start': self.step.masigpro})
        self.masigpro.on('end', self.set_step, {'end': self.step.masigpro})
        self.masigpro.on('end', self.set_output)
        self.masigpro.run()

    def set_output(self):
        shutil.rmtree(self.output_dir)
        shutil.copytree(self.masigpro.output_dir, self.output_dir)
        os.rename(os.path.join(self.output_dir, 'result.tsv'), os.path.join(self.output_dir, 'result.xls'))
        self.option('result').set_path(os.path.join(self.output_dir, 'result.xls'))
        self.set_db()

    def set_db(self):
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/05_Time_series_Analysis')
        self.database = self.api.api('whole_transcriptome.masigpro')
        self.database.add_masigpro(
            result_table=self.option('result').path,
            main_id=self.option('main_id'),
            s3_output=self._sheet.output
        )
        self.end()

    def end(self):
        upload_dir = os.path.join(self.work_dir, 'upload')
        if os.path.isdir(upload_dir):
            os.mkdir(upload_dir)
        os.mkdir(upload_dir)
        for fname in os.listdir(self.output_dir):
            if not fname.endswith('.pkl'):
                os.link(os.path.join(self.output_dir, fname), os.path.join(upload_dir, fname))
        if os.path.exists(os.path.join(upload_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(upload_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(upload_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(upload_dir)
        result_dir.add_relpath_rules([
            [r'.', '', '时序差异分析文件', 0, "211565"]
        ])
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/05 Time_series_Analysis", "", "时序分析结果目录", 0]
        ]

        result_dir.add_regexp_rules([
            [r'result.xls', 'xls', '时序差异分析详情表', 0, "211566"],
            [r'heatmap.*.pdf', 'pdf', '时序差异热图pdf', 0, "211567"],
            [r'heatmap.*.png', 'png', '时序差异热图png', 0, "211568"],
            [r'heatmap.*.svg', 'svg', '时序差异热图svg', 0, "211569"],
            [r'groups.*.pdf', 'pdf', '不同处理在不同时间点的表达模式图pdf', 0, "211570"],
            [r'groups.*.png', 'png', '不同处理在不同时间点的表达模式图png', 0, "211571"],
            [r'groups.*.svg', 'svg', '不同处理在不同时间点的表达模式图svg', 0, "211572"],
            [r'profile.*.pdf', 'pdf', '聚类指定基因在所有样本中的表达模式图pdf', 0, "211573"],
            [r'profile.*.png', 'png', '聚类指定基因在所有样本中的表达模式图png', 0, "211574"],
            [r'profile.*.svg', 'svg', '聚类指定基因在所有样本中的表达模式图svg', 0, "211575"],
            [r'run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(MasigproWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_kmeans(self):
        import random
        from mbio.workflows.ref_rna_v2.report.masigpro import MasigproWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.matrix.tsv',
                'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.design.tsv',
                'cluster': 4,
                'method': 'kmeans'
            }
        }
        wsheet = Sheet(data=data)
        wf = MasigproWorkflow(wsheet)
        wf.sheet.id = 'ref_rna_v2_upgrade'
        wf.sheet.project_sn = 'ref_rna_v2_upgrade'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_kmeans')])
    unittest.TextTestRunner(verbosity=2).run(suite)
