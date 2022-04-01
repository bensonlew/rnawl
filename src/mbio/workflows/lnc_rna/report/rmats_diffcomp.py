# -*- coding: utf-8 -*-
# __author__ = 'shicaiping, qinjincheng'

from biocluster.workflow import Workflow
import os
import shutil
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json


class RmatsDiffcompWorkflow(Workflow):
    '''
    last_modify: 2019.03.19
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsDiffcompWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'delta_psi', 'type': 'float', 'default': None},
            {'name': 'significant_diff', 'type': 'string', 'default': None},
            {'name': 'significant_value', 'type': 'float', 'default': None},
            {'name': 'type_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/05 AS/AS_diff_compare')
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

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        self.get_run_log()
        self.run_download()
        super(RmatsDiffcompWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_splicing_rmats_diffcomp", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_download(self):
        self.step.add_steps('download')
        self.download = self.add_module('lnc_rna.download')
        options = {
            's3_file_list': self.option('s3_file_list'),
            'list_type': 'loc2name'
        }
        self.download.set_options(options)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.run_rmats_diffcomp)
        self.download.run()

    def run_rmats_diffcomp(self):
        self.step.add_steps('rmats_diffcomp')
        self.rmats_diffcomp = self.add_tool('lnc_rna.structure.rmats_diffcomp')
        self.rmats_diffcomp.set_options({
            'rmats_detail_list': self.download.option('local_file_list'),
            'significant_value': self.option('significant_value'),
            'delta_psi': self.option('delta_psi'),
            'significant_diff': self.option('significant_diff'),
            'type_file': self.option('type_file'),
        })
        self.rmats_diffcomp.on('start', self.set_step, {'start': self.step.rmats_diffcomp})
        self.rmats_diffcomp.on('end', self.set_step, {'end': self.step.rmats_diffcomp})
        self.rmats_diffcomp.on('end', self.set_output)
        self.rmats_diffcomp.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.rmats_diffcomp.option('diffcomp_txt').path
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.diffcomp_txt = link_name
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.database = self.api.api('lnc_rna.rmats_diffcomp')
        self.database.add_rmats_diffcomp_detail(
            diffcomp_txt=self.diffcomp_txt,
            main_id=self.option('main_id'),
        )
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.rename(self.output_dir + '/diffcomp.txt', self.output_dir + '/diffcomp_detial.xls')
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["05 AS", "", "可变剪切分析结果目录", 0],
            ["05 AS/AS_diff_compare", "", "差异可变剪切比较文件", 0],
        ]
        relpath = [
            ['.', '', '差异可变剪切事件比较文件', 0],
            ['./diffcomp_detial.xls', 'xls', '差异可变剪切比较结果表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ]
        result_dir.add_relpath_rules(relpath)
        super(RmatsDiffcompWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test(self):
        from mbio.workflows.lnc_rna.report.rmats import RmatsWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'lnc_rna.report.rmats',
            'options': {
                'group_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/group.txt',
                'control_table': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/rmats/control.txt',
                's3_file_list': '/mnt/ilustre/users/sanger-dev/workspace/20190313/SnpIndel_lnc_rna_6850_8188/s3_bam.list',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'read_type': 'paired',
                'lib_type': 'fr-unstranded',
                'control_id': '5c6dfb5517b2bf6465297fca',
                'group_id': '5c6dff2d17b2bf72d19b6db4',
                'submit_location': 'splicingrmats',
                'task_id': 'lnc_rna',
                'task_type': 2,
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsWorkflow(wsheet)
        wf.sheet.id = 'lnc_rna'
        wf.sheet.project_sn = 'lnc_rna'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    unittest.main()