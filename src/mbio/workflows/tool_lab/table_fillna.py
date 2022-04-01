# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""表格合并"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId

class TableFillnaWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TableFillnaWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'whole_transcriptome.common'},  # 上传的表格
            {'name': 'sep', 'type': 'string'},  # 表格的分隔符
            {'name': 'nav', 'type': 'string'},
            {'name': 'rmp_type', 'type': 'bool'},
            {'name': 'rmp', 'type': 'float'},
            {'name': 'method', 'type': 'string'},
            {'name': 'transposed_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        # self.tool = self.add_tool("tool_lab.table_kit_standard")
        self.set_options(self._sheet.options())



    def run(self):
        print self.work_dir
        self.run_fillna()
        super(TableFillnaWorkflow, self).run()

    def run_fillna(self):
        self.fillna = self.add_tool('tool_lab.table_fillna')
        if self.option('rmp_type') == False:
            self.fillna.set_options({
                'table': self.option('table'),
                'sep': self.option('sep'),
                'nav': self.option('nav'),
                'rmp_type': self.option('rmp_type'),
                # 'rmp': self.option('rmp'),
                'method': self.option('method'),
            })
        if self.option('rmp_type') == True:
            self.fillna.set_options({
                'table': self.option('table'),
                'sep': self.option('sep'),
                'nav': self.option('nav'),
                'rmp_type': self.option('rmp_type'),
                'rmp': self.option('rmp'),
                'method': self.option('method'),
            })
        self.fillna.on('end', self.set_output)
        self.fillna.run()

    def set_output(self):
        for file in os.listdir(self.fillna.output_dir):
            os.link(os.path.join(self.fillna.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "excel转化txt文件",0],
        #     ["transposed.txt", "txt", "表格转置文件",0]
        # ])
        super(TableFillnaWorkflow, self).end()
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.table_fillna import TableFillnaWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'fillna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.table_fillna',
            'options': {
                'table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/fillna.txt',
                'sep': 'tab',
                'nav': 'NA',
                'rmp_type': False,
                'method': 'missForest',
                # 'rmp': '0.0'
            }
        }
        wsheet = Sheet(data=data)
        wf =TableFillnaWorkflow(wsheet)
        wf.sheet.id = 'table_fillna'
        wf.sheet.project_sn = 'table_fillna'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
