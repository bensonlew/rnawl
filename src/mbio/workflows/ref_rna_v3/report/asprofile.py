# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
import re
import json
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart
import glob


class AsprofileWorkflow(Workflow):
    """
    可变剪切事件

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AsprofileWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gtf_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'group_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'sample_list', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 AS/05 ASprofile_stat')
        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
                                                  main_id=self.option('main_id'))
            interactiondelete.delete_interactions_records()


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
        super(AsprofileWorkflow, self).send_log(data)


    def check_options(self):
        # for k, v in self.sheet.options().items():
        #     self.logger.debug('{} = {}'.format(k, v))
        pass

    def run(self):
        self.get_run_log()
        self.run_asprofile()
        super(AsprofileWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_asprofile", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_asprofile(self):
        self.asprofile = self.add_module('ref_rna_v3.asprofile')
        self.asprofile.set_options({
            'gtf_dir': self.option('gtf_dir'),
            'ref_fa': self.option('ref_fa'),
            'ref_gtf': self.option('ref_gtf'),
            'group_table': self.option('group_table')
        })
        self.asprofile.on('end', self.set_db)
        self.asprofile.run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        chart.chart_asprofile_stat(os.path.join(self.asprofile.output_dir, 'AS_statistics_merge.txt'))
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            if os.path.exists(self.output_dir + "/" + os.path.basename(p)):
                os.remove(self.output_dir + "/" + os.path.basename(p))
            os.link(p, self.output_dir + "/" + os.path.basename(p))


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.asprofile.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.asprofile.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.asprofile.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.asprofile.output_dir)
        self.inter_dirs = [
            ["04 AS", "", "可变剪切分析结果目录",0],
            ["04 AS/05 ASprofile_stat", "", "可变剪切ASprofile统计结果目录", 0]
        ]

        result_dir.add_relpath_rules([
            [".", "", "可变剪切ASprofile统计文件", 0],
            ['./AS_result_merge.txt', 'txt', 'ASprofile结果详情文件', 0],
            ['./AS_statistics_merge.txt', 'txt', 'ASprofile统计详情文件', 0],
            ['./*.pdf', 'odf', 'ASprofile统计图', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(AsprofileWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        asprofile = self.api.api("ref_rna_v3.asprofile")
        # add result info
        s3_as_result = self._sheet.output + '/' + 'AS_result_merge.txt'
        as_result = os.path.join(self.asprofile.output_dir, 'AS_result_merge.txt')
        as_statistics = os.path.join(self.asprofile.output_dir, 'AS_statistics_merge.txt')
        main_id = asprofile.add_asprofile_result(as_result, s3_as_result, self.option('main_id'))
        asprofile.add_asprofile_statistics(as_statistics, main_id, self.option('sample_list'))
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.asprofile import AsprofileWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "ASprofile" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.asprofile",
            "options": dict(
                gtf_dir='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/gtf/',
                ref_fa='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa'
            )
        }

        wsheet = Sheet(data=data)
        wf =AsprofileWorkflow(wsheet)
        wf.sheet.id = 'ASprofile'
        wf.sheet.project_sn = 'ASprofile'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
