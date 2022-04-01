# -*- coding:utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types


class TransSepWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(TransSepWorkflow, self).__init__(wsheet)
        options = [
                {'name': 'file_name', 'type': 'infile', 'format': 'tool_lab.no_empty'},
            {'name': 'output', 'type': 'string', 'default': 'out.txt'},
            {'name': 'sep_from', 'type': 'string'},
            {'name': 'sep_to', 'type': 'string'},
            {'name': 'main_table', 'type': 'string', 'default': 'trans_sep'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.trans = self.add_tool('tool_lab.trans_sep')

    def run(self):
        self.trans.on('end', self.end)
        self.run_trans()
        super(TransSepWorkflow, self).run()

    def run_trans(self):
        opts = {
            'input': self.option('file_name').path,
            'output': self.option('output'),
            'sep_from': self.option('sep_from'),
            'sep_to': self.option('sep_to'),
        }
        self.trans.set_options(opts)
        self.trans.run()

    def set_db(self):
        remote_dir = self._sheet.output
        file_path = os.path.join(remote_dir, self.option('output'))
        api = self.api.api('tool_lab.api_base')
        query_dict = {'_id': self.option('main_id')}
        update_dict = {'result_path': file_path}
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.trans.output_dir)
        result_dir.add_regexp_rules([
            ['.*.txt', 'txt', '分隔符转换结果文件', 0],
            ['.*.log', 'log', '分隔符转换log日志', 0],
        ])
        result_dir.add_relpath_rules([
            ['.', '', '分隔符转换结果文件夹', 0],
        ])
        super(TransSepWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet

    data = {
            "id": "tt_trans",
            "type":"workflow",
            "name": "test_trans_sep",
            "options":{
                "file_name": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/test/t.txt",
                "output": "out",
                "sep_from": "semicolon",
                "sep_to": "tab"
                }
            
            }
    wsheet = Sheet(data=data)
    wf = TransSepWorkflow(wsheet)
    wf.run()

