# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import time
import pandas as pd

class DeletedExpbatchWorkflow(Workflow):
    """
    可变剪切事件

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DeletedExpbatchWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'library', 'type': 'string'},
            {'name': 'level', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))


    def run(self):
        self.start_listener()
        self.fire("start")
        self.set_db()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        time.sleep(7)

        deleted_expbatch = self.api.api("whole_transcriptome.deleted_expbatch")
        deleted_expbatch.deleted_batch(task_id=self.option('task_id'), library=self.option('library'), level=self.option('level'), main_id=self.option('main_id'))
        self.end()



    def end(self):
        super(DeletedExpbatchWorkflow, self).end()

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
            "id": "add_annotation" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "ref_rna_v3.add_annotation",
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
