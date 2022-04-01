# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.whole_transcriptome.functions import workfuncdeco
import os
import shutil
import unittest
import re
import json
import pandas as pd
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.whole_transcriptome.chart.chart import Chart
import glob
from biocluster.config import Config
from mbio.packages.project_demo.interaction_rerun.interaction_delete import linkfile

class RmatsStatWorkflow(Workflow):
    '''
    last_modify: 2019.06.14
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsStatWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'root', 'type': 'infile', 'format': 'whole_transcriptome.common_dir'},
            {'name': 'pvalue_fdr', 'type': 'string', 'default': None},
            {'name': 'fdr', 'type': 'float', 'default': None},
            {'name': 'psi', 'type': 'float', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
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
        super(RmatsStatWorkflow, self).send_log(data)

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @workfuncdeco
    def run(self):
        self.get_run_log()
        self.run_rmats_stat()
        super(RmatsStatWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="splicing_rmats_stats", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_rmats_stat(self):
        self.step.add_steps('rmats_stat')
        self.rmats_stat = self.add_tool('whole_transcriptome.structure.rmats_stat')
        self.rmats_stat.set_options({
            'root': self.option('root'),
            'pvalue_fdr': self.option('pvalue_fdr'),
            'fdr': self.option('fdr'),
            'psi': self.option('psi')
        })
        self.rmats_stat.on('start', self.set_step, {'start': self.step.rmats_stat})
        self.rmats_stat.on('end', self.set_step, {'end': self.step.rmats_stat})
        self.rmats_stat.on('end', self.set_output)
        self.rmats_stat.run()

    def set_output(self, event):
        for src in [self.rmats_stat.option('event_stats').path, self.rmats_stat.option('psi_stats').path]:
            dst = os.path.join(self.output_dir, os.path.basename(src))
            shutil.copy(src, dst)
            self.logger.info('succeed in copying {} to {}'.format(src, dst))
        else:
            self.set_db()

    @workfuncdeco
    def set_db(self):
        api = self.api.api('whole_transcriptome.rmats')
        api.add_splicing_rmats_stats_for_controller(stat_id=self.option('main_id'), outpath=self.output_dir)
        self.end()

    def chart(self):
        for (key, value) in [["LD_LIBRARY_PATH",Config().SOFTWARE_DIR + "/bioinfo/sg_chart/miniconda2/lib:$LD_LIBRARY_PATH"],["NODE_PATH",Config().SOFTWARE_DIR + "/bioinfo/sg_chart/node-v14.16.0-linux-x64/lib/node_modules"]]:
            if key not in os.environ.keys():
                os.environ[key] = value
            else:
                os.environ[key] = value + ":" + os.environ[key]
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        splice_diff = self.output_dir + "/event_stats.file.txt"
        splice_psi = self.output_dir + "/psi_stats.file.txt"
        chart.chart_splice_diff_stat(splice_diff, splice_psi, cmp_name=self.option("root").prop["path"].split("/")[-2])
        chart.to_pdf()

        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            linkfile(p, self.output_dir + "/" + os.path.basename(p))

    @workfuncdeco
    def end(self):
        self.chart()
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/05 AS/AS_diff_stat')
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 AS", "", "可变剪切分析结果目录",0],
            ["05 AS/AS_diff_stat", "", "差异可变剪切统计文件", 0]
        ]
        relpath = [
            ['.', '', '差异可变剪切事件统计',0],
            ['event_stats.file.txt', '', '差异可变剪切事件统计表',0],
            ['psi_stats.file.txt', '', '差异可变剪切模式变化统计表',0],
            ['*.pdf', 'txt', "差异可变剪切统计图", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsStatWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test(self):
        from mbio.workflows.whole_transcriptome.report.rmats_stat import RmatsStatWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.report.rmats_stat',
            'options': {
                'root': '/mnt/ilustre/users/sanger-dev/workspace/20190719/Refrna_tsg_34902/upload/09AS/P4_MG_vs_P4_Br3',
                'pvalue_fdr': 'fdr',
                'fdr': '0.06',
                'psi': '0.1',
                'main_id': '5d90139517b2bf3331919f6f'
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsStatWorkflow(wsheet)
        wf.sheet.id = 'whole_transcriptome'
        wf.sheet.project_sn = 'whole_transcriptome'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
