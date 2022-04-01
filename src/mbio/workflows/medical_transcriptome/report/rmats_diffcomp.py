# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.medical_transcriptome.functions import workfuncdeco
import os
import shutil
import unittest
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart


class RmatsDiffcompWorkflow(Workflow):
    '''
    last_modify: 2019.06.14
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsDiffcompWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'delta_psi', 'type': 'float', 'default': None},
            {'name': 'significant_diff', 'type': 'string', 'default': None},
            {'name': 'significant_value', 'type': 'float', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 Gene_structure_analysis/01 AS/03 AS_diff_compare')
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
        super(RmatsDiffcompWorkflow, self).send_log(data)

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
        self.run_download()
        super(RmatsDiffcompWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_splicing_rmats_diffcomp", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_download(self):
        self.step.add_steps('download')
        self.download = self.add_module('medical_transcriptome.download')
        self.download.set_options({'s3_file_list': self.option('s3_file_list')})
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.run_rmats_diffcomp)
        self.download.run()

    def run_rmats_diffcomp(self):
        self.step.add_steps('rmats_diffcomp')
        self.rmats_diffcomp = self.add_tool('medical_transcriptome.structure.rmats_diffcomp')
        self.rmats_diffcomp.set_options({
            'rmats_detail_list': self.download.option('my_file_list'),
            'delta_psi': self.option('delta_psi'),
            'significant_diff': self.option('significant_diff'),
            'significant_value': self.option('significant_value')
        })
        self.rmats_diffcomp.on('start', self.set_step, {'start': self.step.rmats_diffcomp})
        self.rmats_diffcomp.on('end', self.set_step, {'end': self.step.rmats_diffcomp})
        self.rmats_diffcomp.on('end', self.set_output)
        self.rmats_diffcomp.run()

    def set_output(self, event):
        source = self.rmats_diffcomp.option('diffcomp_txt').path
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.diffcomp_txt = link_name
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.set_db()

    @workfuncdeco
    def set_db(self):
        api = self.api.api('medical_transcriptome.rmats_diffcomp')
        api.add_rmats_diffcomp_detail(diffcomp_txt=self.diffcomp_txt, main_id=self.option('main_id'))
        self.end()



    @workfuncdeco
    def end(self):
        # if os.path.exists(os.path.join(self.rmats_diffcomp.output_dir, os.path.basename(self.run_log))):
        #     os.remove(os.path.join(self.rmats_diffcomp.output_dir, os.path.basename(self.run_log)))
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/01 AS", "", "可变剪切分析结果目录", 0],
            ["03 Gene_structure_analysis/01 AS/03 AS_diff_compare", "", "差异可变剪切比较文件", 0]
        ]
        relpath = [
            ['.', '', '差异可变剪切事件比较文件',0,"211641"],
            ['./diffcomp.txt', 'txt', '差异可变剪切比较结果表', 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsDiffcompWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test(self):
        from mbio.workflows.medical_transcriptome.report.rmats_diffcomp import RmatsDiffcompWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'medical_transcriptome.report.rmats_diffcomp',
            'options': {
                's3_file_list': '/mnt/ilustre/users/sanger-dev/workspace/20200812/RmatsDiffcomp_tsg_38314_5063_7579/table.list',
                'delta_psi': '0.0',
                'significant_diff': 'fdr',
                'significant_value': '0.05',
                'main_id': '5f33b28a17b2bf7043eec5e0'
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsDiffcompWorkflow(wsheet)
        wf.sheet.id = 'tsg_33555_8908_5404'
        wf.sheet.project_sn = '188_5c820f0e8b599'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
