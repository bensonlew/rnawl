# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import json
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class DiffAsprofileWorkflow(Workflow):
    """
    可变剪切事件

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffAsprofileWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "AS_result_merge", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': 'sample', 'type': 'string'},
            {'name': 'group_dict', 'type': 'string'},
            {'name': 'filter', 'type': 'int'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/03 Gene_structure_analysis/01 AS/06 ASprofile_diff')
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
        super(DiffAsprofileWorkflow, self).send_log(data)

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))

    def run(self):
        self.get_run_log()
        self.run_diff_as()
        super(DiffAsprofileWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_asprofile_diff", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_diff_as(self):
        self.diff_as = self.add_tool('medical_transcriptome.ASprofile.diff_asprofile')
        if self.option('sample'):
            self.diff_as.set_options({
                'AS_result_merge': self.option('AS_result_merge'),
                'sample': self.option('sample'),
                'filter': self.option('filter')
            })
        if self.option('group_dict'):
            self.diff_as.set_options({
                'AS_result_merge': self.option('AS_result_merge'),
                'group_dict': self.option('group_dict'),
                'filter': self.option('filter')
            })
        self.diff_as.on('end', self.set_db)
        self.diff_as.run()

    def end(self):
        if os.path.exists(os.path.join(self.diff_as.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.diff_as.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.diff_as.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.diff_as.output_dir)
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/01 AS", "", "可变剪切分析结果目录", 0],
            ["03 Gene_structure_analysis/01 AS/06 ASprofile_diff", "", "可变剪切ASprofile差异结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "可变剪切ASprofile差异文件", 0],
            ['./AS_diff.txt', 'txt', 'ASprofile差异结果详情文件'],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(DiffAsprofileWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        asprofile = self.api.api("medical_transcriptome.asprofile")
        # add result info
        diff_asprofile = os.path.join(self.diff_as.output_dir, 'AS_diff.txt')
        asprofile.add_diff_asprofile_result(diff_asprofile, self.option('main_id'))
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.diff_asprofile import DiffAsprofileWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "Diff_ASprofile" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.report.diff_asprofile",
            "options": dict(
                AS_result_merge='/mnt/ilustre/users/sanger-dev/workspace/20200813/Asprofile_tsg_38314_8419_7589/Asprofile/output/AS_result_merge.txt',
                sample='H1581_1,H1581_2,H1581_3',
                filter=2
            )
        }

        wsheet = Sheet(data=data)
        wf =DiffAsprofileWorkflow(wsheet)
        wf.sheet.id = 'diff_asprofile'
        wf.sheet.project_sn = 'diff_asprofile'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
