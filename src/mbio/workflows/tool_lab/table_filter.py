# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""表格筛选"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId

class TableFilterWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TableFilterWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'sep', 'type': 'string'},
            {'name': 'table_head', 'type': 'bool'},
            {'name': 'contain', 'type': 'bool'},
            {'name': 'lower', 'type': 'bool'},
            {'name': 'number', 'type': 'int'},
            {'name': 'info', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'filter_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        # self.tool = self.add_tool("tool_lab.table_standard")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_filter()
        super(TableFilterWorkflow, self).run()

    def run_filter(self):
        self.table_filter = self.add_tool('tool_lab.table_filter')
        self.table_filter.set_options({
            'table': self.option('table'),
            'sep': self.option('sep'),
            'table_head': self.option('table_head'),
            'contain': self.option('contain'),
            'lower': self.option('lower'),
            'number': self.option('number'),
            'info': self.option('info'),
            'filter_table': self.option('filter_table')
        })
        self.table_filter.on('end', self.set_output)
        self.table_filter.run()

    def set_output(self):
        p = self.table_filter.option('filter_table').path

        link_names = os.path.join(self.output_dir, os.path.basename(p))
        os.link(p, link_names)
        self.option('filter_table').set_path(link_names)
        self.end()
    # def set_db(self):
    #     """
    #     保存结果标准化数据到mongo数据库中
    #     """
    #     # record_id = self.option("standard_id")
    #     # if isinstance(record_id, types.StringTypes):
    #     #     record_id = ObjectId(record_id)
    #     # elif isinstance(record_id, ObjectId):
    #     #     record_id = record_id
    #     # else:
    #     #     self.set_error("main_id参数必须为字符串或者ObjectId类型!", code="13701601")
    #     standard = self.api.api("tool_lab.table_kit")
    #     # add result info
    #     print self.option('transposed_table').path
    #     standard.add_table_standard()
    #     self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "标准化分析文件",0],
            ["standard.txt", "txt", "标准化分析",0]
        ])
        super(TableFilterWorkflow, self).end()
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.table_standard import TableStandardWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'standard_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.table_standard',
            'options': {
                'sep' : 'tab',
                'observe' : 'sum',
                'feature' : 'standard_scale',
                'table' : '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/table.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf =TableStandardWorkflow(wsheet)
        wf.sheet.id = 'table_standard'
        wf.sheet.project_sn = 'table_standard'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
