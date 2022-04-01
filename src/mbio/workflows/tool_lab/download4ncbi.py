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


class Download4ncbiWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Download4ncbiWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.run_download()
        super(Download4ncbiWorkflow, self).run()

    def run_download(self):
        self.download = self.add_tool("tool_lab.download4ncbi.download")
        self.download.set_options({
            'id': self.option('id'),
        })
        self.download.on('end', self.set_db)
        self.download.run()

    def set_db(self):

        """
         保存结果标准化数据到mongo数据库中
        """
        download = self.api.api('tool_lab.download4ncbi')
        download.add_download(main_id=self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.download.output_dir)
        super(Download4ncbiWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.download4ncbi import Download4ncbiWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'downncbi_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.download4ncbi',
            'options': {
                'id': 'GCF_004214875.1',
            }
        }
        wsheet = Sheet(data=data)
        wf =Download4ncbiWorkflow(wsheet)
        wf.sheet.id = 'download'
        wf.sheet.project_sn = 'download'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
