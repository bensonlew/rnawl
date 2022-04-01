# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class DiffStringPicturesWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffStringPicturesWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'string_xml', 'type': 'infile', 'format': 'itraq_and_tmt.common'},
            {'name': 'gene_list', 'type': 'string'},
            {'name': 'identity', 'type': 'float', 'default': 98},
            {'name': 'max_num', 'type': 'int', 'default': 300},
            {'name': 'species', 'type': 'int', 'default': 0},
            {'name': 'sepcies_classification', 'type': 'string'},
            {'name': 'useblast', 'type': 'string', 'default': 'no'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def run(self):
        self.run_string()
        super(DiffStringPicturesWorkflow, self).run()

    def run_string(self):
        self.string = self.add_tool("tool_lab.diff_string_pictures")
        self.string.set_options({
            'gene_list': self.option('gene_list'),
            'string_xml': self.option('string_xml'),
            'identity': self.option('identity'),
            'max_num': self.option('max_num'),
            'species': self.option('species'),
            'useblast': self.option('useblast')
        })
        self.string.on('end', self.set_db)
        self.string.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        s3_output = os.path.join(self._sheet.output, 'protein.list/protein.list.network.png')
        string_pictures = self.api.api('tool_lab.string_pictures')
        output_dir = os.path.join(self.string.output_dir, 'protein.list')
        string_pictures.add_string(main_id=self.option('main_id'), string_dir=output_dir, s3_path=s3_output)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.string.output_dir)
        super(DiffStringPicturesWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.batch import BatchWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'batch_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.batch',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/unigene.tpm.matrix.annot.txt',
                'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/group_table.txt',
                'batch_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/batch/batch_new.txt',
                'batch_method': 'combat',
                'has_batch': True
            }
        }
        wsheet = Sheet(data=data)
        wf =BatchWorkflow(wsheet)
        wf.sheet.id = 'batch_effect'
        wf.sheet.project_sn = 'batch_effect'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
