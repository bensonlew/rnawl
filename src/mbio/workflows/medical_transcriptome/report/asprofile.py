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
from mbio.packages.medical_transcriptome.chart.chart import Chart
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
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 Gene_structure_analysis/01 AS/05 ASprofile_stat')
        self.logger.info("klalalala{}".format(self._sheet.output))
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
        get_run_log = GetRunLog("medical_transcriptome", table="sg_asprofile", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_asprofile(self):
        self.asprofile = self.add_module('medical_transcriptome.asprofile.asprofile')
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
            os.link(p, self.asprofile.output_dir + "/" + os.path.basename(p))

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.asprofile.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.asprofile.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.asprofile.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.asprofile.output_dir)
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/01 AS", "", "可变剪切分析结果目录", 0],
            ["03 Gene_structure_analysis/01 AS/05 ASprofile_stat", "", "可变剪切ASprofile统计结果目录", 0]
        ]

        result_dir.add_relpath_rules([
            [".", "", "可变剪切ASprofile统计文件", 0],
            ['./AS_result_merge.txt', 'txt', 'ASprofile结果详情文件', 0],
            ['./AS_statistics_merge.txt', 'txt', 'ASprofile统计详情文件', 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*\.column\.pdf', 'pdf', 'ASprofile统计堆积柱形图', 0],
            [r".*.\pie\.pdf", 'pdf', "ASprofile单样本统计饼图", 0]
        ])
        super(AsprofileWorkflow, self).end()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        asprofile = self.api.api("medical_transcriptome.asprofile")
        # add result info
        s3_as_result = self._sheet.output + '/' + 'AS_result_merge.txt'
        as_result = os.path.join(self.asprofile.output_dir, 'AS_result_merge.txt')
        as_statistics = os.path.join(self.asprofile.output_dir, 'AS_statistics_merge.txt')
        main_id = asprofile.add_asprofile_result(as_statistics, s3_as_result, self.option('main_id'))

        asprofile.add_asprofile_statistics(as_statistics, main_id, self.option('sample_list'))
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.asprofile import AsprofileWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "ASprofile" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.asprofile",
            "output":"s3://medical_transcriptome/files/test/medical_transcriptome/medical_transcriptome/interaction_results/ASprofile_20200813_092731",
            "options": dict(
                gtf_dir='s3://refrnav2/files/m_188/188_5d01dede4f911/tsg_38314/intermediate_results/Stringtie/',
                ref_fa='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                ref_gtf = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/gtf/Homo_sapiens.GRCh38.96.gtf",
                group_table="/mnt/ilustre/users/sanger-dev/workspace/20200813/Asprofile_tsg_38314_8419_7589/group_table.txt",
                sample_list="H1581_1,H1581_2,H1581_3,H1581_4,H1581_5,H1581_6,H1581_7,H1581_8,H1581_9,SNU16_1,SNU16_2,SNU16_3,SNU16_4,SNU16_5,SNU16_6,SNU16_7,SNU16_8,SNU16_9",
                main_id="5f34970317b2bf7920d32455"
            )
        }

        wsheet = Sheet(data=data)
        wf =AsprofileWorkflow(wsheet)
        wf.sheet.id = 'medical_transcriptome'
        wf.sheet.project_sn = 'medical_transcriptome'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
