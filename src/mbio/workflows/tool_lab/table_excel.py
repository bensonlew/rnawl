# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""表格合并"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId

class TableExcelWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TableExcelWorkflow, self).__init__(wsheet_object)
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
        self.run_excel()
        super(TableExcelWorkflow, self).run()

    def run_excel(self):
        self.excel = self.add_tool('tool_lab.table_excel')
        self.excel.set_options({
            'table': self.option('table'),
            'sep': self.option('sep'),
        })
        self.excel.on('end', self.set_output)
        self.excel.run()

    def set_output(self):
        for file in os.listdir(self.excel.output_dir):
            os.link(os.path.join(self.excel.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "excel转化txt文件",0],
        #     ["transposed.txt", "txt", "表格转置文件",0]
        # ])
        super(TableExcelWorkflow, self).end()
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.table_excel import TableExcelWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'excel_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.table_excel',
            'options': {
                'sep' : 'tab',
                'table' : '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/test_excel.xlsx',
            }
        }
        wsheet = Sheet(data=data)
        wf =TableExcelWorkflow(wsheet)
        wf.sheet.id = 'table_excel'
        wf.sheet.project_sn = 'table_excel'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
