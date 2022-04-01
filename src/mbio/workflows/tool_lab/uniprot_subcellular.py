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


class UniprotSubcellularWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(UniprotSubcellularWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'uniprot_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.run_uniprot_api()
        super(UniprotSubcellularWorkflow, self).run()

    def run_uniprot_api(self):
        self.uniprot_api = self.add_tool("tool_lab.uniprot_subcellular")
        self.uniprot_api.set_options({
            'uniprot_id': self.option('uniprot_id')
        })
        self.uniprot_api.on('end', self.set_db)
        self.uniprot_api.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.uniprot_api.output_dir)
        super(UniprotSubcellularWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.protein_picture import ProteinPictureWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'picture_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.protein_picture',
            'options': {
                'picture': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/protein_picture/Dl09YYJ_bjhb1_1580_Result.PNG',
                'if_black': 'no',
            }
        }
        wsheet = Sheet(data=data)
        wf =ProteinPictureWorkflow(wsheet)
        wf.sheet.id = 'protein_picture'
        wf.sheet.project_sn = 'protein_picture'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
