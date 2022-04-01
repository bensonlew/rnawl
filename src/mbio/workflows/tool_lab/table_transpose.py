# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""表格合并"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId

class TableTransposeWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TableTransposeWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'whole_transcriptome.common'},  # 上传的表格
            {'name': 'sep', 'type': 'string'},  # 表格的分隔符
            {'name': 'transposed_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        # self.tool = self.add_tool("tool_lab.table_kit_standard")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_transpose()
        super(TableTransposeWorkflow, self).run()

    def run_transpose(self):
        self.transpose = self.add_tool('tool_lab.table_transpose')
        self.transpose.set_options({
            'table': self.option('table'),
            'sep': self.option('sep'),
        })
        self.transpose.on('end', self.set_output)
        self.transpose.run()

    def set_output(self):
        p = self.transpose.option('transposed_table').path

        link_names = os.path.join(self.output_dir, os.path.basename(p))
        os.link(p, link_names)
        self.option('transposed_table').set_path(link_names)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表格转置文件",0],
            ["transposed.txt", "txt", "表格转置文件",0]
        ])
        super(TableTransposeWorkflow, self).end()
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.table_transpose import TableTransposeWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'transpose_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.table_transpose',
            'options': {
                'sep' : 'tab',
                'table' : '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/venn.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf =TableTransposeWorkflow(wsheet)
        wf.sheet.id = 'table_transpose'
        wf.sheet.project_sn = 'table_transpose'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
