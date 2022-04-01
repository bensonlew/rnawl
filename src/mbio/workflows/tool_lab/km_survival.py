# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
import glob
import os
import json
import time
import re
from biocluster.file import getsize, exists
from biocluster.file import download
import pandas as pd
import unittest


class KmSurvivalWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(KmSurvivalWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'km_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'risk_table', 'type': 'bool', 'default': True},
            {'name': 'conf', 'type': 'bool', 'default': True},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("tool_lab.km_survival")
        self.km = self.api.api("tool_lab.km")

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(KmSurvivalWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        pdf_path = self.tool.option('km_pdf').path
        s3_output = '{}/{}'.format(self._sheet.output, os.path.basename(pdf_path))

        # add result info
        self.km.add_km(
                                 main_id=self.option('main_id'),
                                 s3_output=s3_output,
                                 )
        self.set_output()
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "差异分析结果目录", 0, "211063"],
        # ])
        # result_dir.add_regexp_rules([
        #     [r'.*_vs_.*\.xls', 'xls', '差异表达基因详情表',0,"211518"],
        #     [r'.*summary.*\.xls', 'xls', '差异表达基因统计表',0,"211519"],
        #     [r'.*total_diff_stat.*\.xls', 'xls', '差异表达基因详情总表',0,"211520"],
        # ])
        super(KmSurvivalWorkflow, self).end()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
    def run_tool(self):
        opts = dict(
            km_table=self.option('km_table'),
            risk_table=self.option('risk_table'),
            conf=self.option('conf')
        )
        self.tool.set_options(opts)
        self.tool.run()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.diffexp_deseq2 import DiffexpDeseq2Workflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.diffexp_deseq2',
            "options": dict(
                count='/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37857/Quant/output/gene.count.matrix',
                group="/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37857/remote_input/group_table/group.txt",
                cmp= "/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37857/remote_input/control_file/compare.txt",
                padjust_way="BH",
                pvalue=0.05,
                fc=2,
                test='Wald'


            )
        }
        wsheet = Sheet(data=data)
        wf = DiffexpDeseq2Workflow(wsheet)
        wf.sheet.project_sn = '188_5dba6f542345b'
        wf.sheet.task_id = '123456'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


