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

class AddAnnotationWorkflow(Workflow):
    """
    可变剪切事件

    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AddAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'change_database', 'type': 'string'},
            {'name': 'annotation_file', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'annotation_type', 'type': 'string'}, #[replace, supplement]
            {'name': 'annotation_method', 'type': 'string', 'default': 'all'},#[all, name, description]
            {'name': 'update_info', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
            {'name': 'ghost', 'type': 'bool', 'default': False}
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

        add_annotation = self.api.api("ref_rna_v3.add_annotation")
        if self.option('ghost') is True:
            add_annotation.ghost_add_annotation(main_id=self.option('main_id'), task_id=self.option('task_id'))
        else:
            if self.option('change_database') == 'custom':
                annotation_path = os.path.join(self.work_dir, 'annotation_file.txt')
                annotation_file = pd.read_table(self.option('annotation_file').prop['path'], header=None, sep='\t')
                annotation_file.drop_duplicates(subset=[0], keep='first', inplace=True)
                annotation_file.to_csv(annotation_path, index=False, header=False,sep='\t')
                add_annotation.add_name_description(annotation_path, self.option('annotation_type'),
                                                    self.option('annotation_method'), main_id=self.option('main_id'),
                                                    task_id=self.option('task_id'))
            else:
                add_annotation.add_name_description_database(self.option('change_database').lower(), self.option('annotation_type'),
                                                             self.option('annotation_method'), main_id=self.option('main_id'),
                                                             task_id=self.option('task_id'))
        self.end()



    def end(self):
        super(AddAnnotationWorkflow, self).end()

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
