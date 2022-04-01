# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

"""表格合并"""

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
import zipfile

class TableCombineWorkflow(Workflow):
    """
    处理表格
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TableCombineWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'table1', 'type': 'infile', 'format': 'whole_transcriptome.common'},  # 上传的表格
            {'name': 'table_zip', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'sep', 'type': 'string'},  # 表格的分隔符
            {'name': 'method', 'type': 'string'},
            {'name': 'fill', 'type': 'string'},
            {'name': 'fill_value', 'type': 'string'},
            {'name': 'transposed_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        # self.tool = self.add_tool("tool_lab.table_kit_standard")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_combine()
        super(TableCombineWorkflow, self).run()

    def run_unzip(self):
        table_list = list()
        table_list.append(self.option('table1').path)
        filename = os.path.basename(self.option('table_zip').path)
        if filename[-4:] == '.zip':
            fz = zipfile.ZipFile(self.option('table_zip').path, 'r')
            fz_path = os.path.join(self.work_dir, 'zip')
            for f in fz.namelist():
                fz.extract(f, fz_path)

            for root, dirs, files in os.walk(fz_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    table_list.append(file_path)
        else:
            table_list.append(self.option('table_zip').path)
        table_list_str = ','.join(table_list)
        return table_list_str
    def run_combine(self):
        self.combine = self.add_tool('tool_lab.table_combine')
        table_list = self.run_unzip()
        if self.option('fill') == 'other':
            self.combine.set_options({
                'table_list': table_list,
                'sep': self.option('sep'),
                'method':self.option('method'),
                'fill':self.option('fill_value'),
            })
        else:
            if self.option('fill') == '--':

                self.combine.set_options({
                    'table_list': table_list,
                    'sep': self.option('sep'),
                    'method':self.option('method'),
                    'fill':'doubleline',
                })
            else:
                self.combine.set_options({
                    'table_list': table_list,
                    'sep': self.option('sep'),
                    'method':self.option('method'),
                    'fill':self.option('fill'),
                })
        self.combine.on('end', self.set_output)
        self.combine.run()

    def set_output(self):
        p = self.combine.option('transposed_table').path

        link_names = os.path.join(self.output_dir, os.path.basename(p))
        os.link(p, link_names)
        self.option('transposed_table').set_path(link_names)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表格合并文件",0],
            ["combine.txt", "txt", "表格合并文件",0]
        ])
        super(TableCombineWorkflow, self).end()
class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.table_combine import TableCombineWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'combine_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.table_combine',
            'options': {
                'sep' : 'tab',
                'method' : 'all',
                'fill' : 'other',
                'fill_value': 'yes',
                'table1' : '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/combine.txt',
                'table_zip': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/combine1.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf =TableCombineWorkflow(wsheet)
        wf.sheet.id = 'table_combine'
        wf.sheet.project_sn = 'table_combine'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
