# -*- coding: utf-8 -*-
# __author__ = 'shicaiping,qinjincheng'

from biocluster.workflow import Workflow
from mbio.packages.whole_transcriptome.functions import workfuncdeco
import os
import shutil
from mbio.packages.whole_transcriptome.upload import Upload
# from mbio.packages.whole_transcriptome.upload import get_gid2des_dct, export_rmats_detail, export_rmats_diff_stats, export_rmats_psi
import re
import json
import pandas as pd
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class RmatsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RmatsWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 's3_file_list', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'group_table', 'type': 'infile', 'format': 'sample.group_table'},
            {'name': 'control_table', 'type': 'infile', 'format': 'sample.control_table'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'whole_transcriptome.gtf'},
            # ['paired', 'single']
            {'name': 'seq_type', 'type': 'string', 'default': 'paired'},
            # ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']
            {'name': 'lib_type', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'splicing_id', 'type': 'string', 'default': None},
            {'name': 'task_id', 'type': 'string', 'default': None},
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
        super(RmatsWorkflow, self).send_log(data)

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
        super(RmatsWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="splicing_rmats", main_id=self.option('splicing_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_download(self):
        self.step.add_steps('download')
        self.download = self.add_module('whole_transcriptome.download')
        options = {
            's3_file_list': self.option('s3_file_list')
        }
        self.download.set_options(options)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.run_rmats)
        self.download.run()

    def run_rmats(self):
        self.step.add_steps('rmats')
        self.rmats = self.add_module('whole_transcriptome.rmats')
        self.rmats.set_options({
            'control_table': self.option('control_table'),
            'group_table': self.option('group_table'),
            'bam_input': self.download.output_dir,
            'input_gtf': self.option('ref_gtf'),
            'seq_type': self.option('seq_type'),
            'lib_type': self.option('lib_type')
        })
        self.rmats.on('start', self.set_step, {'start': self.step.rmats})
        self.rmats.on('end', self.set_step, {'end': self.step.rmats})
        self.rmats.on('end', self.set_output)
        self.rmats.run()

    def set_output(self, event):
        shutil.rmtree(self.output_dir)
        shutil.copytree(self.rmats.output_dir, self.output_dir)
        self.set_db()

    @workfuncdeco
    def set_db(self):
        api = self.api.api('whole_transcriptome.rmats')
        num, vs_list = self.option('control_table').get_control_info()
        group_spname = self.option('group_table').get_group_spname()
        ctrl, test = vs_list[0]
        self.stat_id = api.add_splicing_rmats_for_controller(
            splicing_id=self.option('splicing_id'),
            outpath=os.path.join(self.output_dir, '{}_vs_{}'.format(ctrl, test)),
            s3_output=self._sheet.output
        )
        self.logger.info("stat is is{}".format(self.stat_id))
        self.set_upload()

    @workfuncdeco
    def set_upload(self):
        self.target_dir = os.path.join(self.work_dir, 'upload')
        if os.path.isdir(self.target_dir):
            shutil.rmtree(self.target_dir)
        os.mkdir(self.target_dir)
        set_upload = Upload()
        if self.option("task_id"):
            gid2des_dct,gid2name_dct = set_upload.get_gid2des_dct(self.option("task_id"))
        else:
            gid2des_dct,gid2name_dct = set_upload.get_gid2des_dct(self.sheet.task_id)
        num, vs_list = self.option('control_table').get_control_info()
        s1, s2 = vs_list[0]
        output_dir = os.path.join(self.target_dir, '{}_vs_{}'.format(s1, s2))
        set_upload.export_rmats_detail(self.option('splicing_id'), gid2des_dct,gid2name_dct, s2, s1, output_dir)
        set_upload.export_rmats_diff_stats(self.stat_id, output_dir)
        set_upload.export_rmats_psi(self.stat_id, output_dir)
        self.end()

    @workfuncdeco
    def end(self):
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/05 AS/AS_diff')
        # self.upload_to_s3(self.output_dir, self._sheet.output)
        if os.path.exists(os.path.join(self.target_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.target_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.target_dir, os.path.basename(self.run_log)))
        if not self.target_dir.endswith("/"):
            self.target_dir = self.target_dir + "/"
        result_dir = self.add_upload_dir(self.target_dir)
        self.inter_dirs = [
            ["05 AS", "", "可变剪切分析结果目录",0],
            ["05 AS/AS_diff", "", "差异可变剪切分析文件", 0]
        ]
        result_dir.add_relpath_rules([
            [r'.', '', '差异可变剪切分析结果文件',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*', '', '差异组别',0,"211627"],
            [r'.*_vs_.*/JC/SE.detail.xls', '', 'SE可变剪切时间详情表（JC）',0],
            [r'.*_vs_.*/JCEC/SE.detail.xls', '', 'SE可变剪切时间详情表（JCEC）',0],
            [r'.*_vs_.*/JC/MXE.detail.xls', '', 'MXE可变剪切时间详情表（JC）',0],
            [r'.*_vs_.*/JCEC/MXE.detail.xls', '', 'MXE可变剪切时间详情表（JCEC）',0],
            [r'.*_vs_.*/JC/A3SS.detail.xls', '', 'A3SS可变剪切时间详情表（JC）',0],
            [r'.*_vs_.*/JCEC/A3SS.detail.xls', '', 'A3SS可变剪切时间详情表（JCEC）',0],
            [r'.*_vs_.*/JC/A5SS.detail.xls', '', 'A5SS可变剪切时间详情表（JC）',0],
            [r'.*_vs_.*/JCEC/A5SS.detail.xls', '', 'A5SS可变剪切时间详情表（JCEC）',0],
            [r'.*_vs_.*/JC/RI.detail.xls', '', 'RI可变剪切时间详情表（JC）',0],
            [r'.*_vs_.*/JCEC/RI.detail.xls', '', 'RI可变剪切时间详情表（JCEC）',0],
            [r'.*_vs_.*/diff_event_stats.xls', '', '组内差异可变剪切事件统计表',0],
            [r'.*_vs_.*/diff_pattern_stats.JC.xls', '', '组内差异可变剪切模式变化统计表（JC）',0],
            [r'.*_vs_.*/diff_pattern_stats.JCEC.xls', '', '组内差异可变剪切模式变化统计表（JCEC）',0]
        ])
        super(RmatsWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''
    def test(self):
        from mbio.workflows.ref_rna_v3.report.rmats import RmatsWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v3.report.rmats',
            'options': {
                's3_file_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/wftest/s3.bam.list',
                'control_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/wftest/control.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/wftest/group.txt',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/rmats/ref_and_new.gtf',
                'seq_type': 'paired',
                'lib_type': 'fr-unstranded'
            }
        }
        wsheet = Sheet(data=data)
        wf = RmatsWorkflow(wsheet)
        wf.sheet.id = 'ref_rna_v3'
        wf.sheet.project_sn = 'ref_rna_v3'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
