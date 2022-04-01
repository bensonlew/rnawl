# -*- coding: utf-8 -*-
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


class FastaFromIdWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(FastaFromIdWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'fasta_file', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'id_file', 'type': 'infile', 'format': 'tool_lab.no_empty'},
            {'name': 'main_table', 'type': 'string', 'default': 'fasta_from_id'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.fasta = self.add_tool('tool_lab.fasta_from_id')

    def run(self):
        self.fasta.on('end', self.end)
        self.run_fasta()
        super(FastaFromIdWorkflow, self).run()

    def run_fasta(self):
        opts = {
            'in_seq': self.option('fasta_file').prop['path'],
            'id_list': self.option('id_file').path,
            'out_seq': 'out.fasta',
        }
        self.fasta.set_options(opts)
        self.fasta.run()

    def set_db(self):
        #remote_dir = self._sheet.output
        #file_path = os.path.join(remote_dir, self.option('output'))
        #api = self.api.api('tool_lab.api_base')
        #query_dict = {'_id': self.option('main_id')}
        #update_dict = {'result_path': file_path}
        #api.update_db_record(self.option('main_table'), query_dict, update_dict)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.fasta.output_dir)
        result_dir.add_regexp_rules([
            [r'.*.fasta', 'fasta', 'fasta提取结果', 0],
            [r'.*.fa', 'fa', 'fasta提取结果', 0],
        ])
        result_dir.add_relpath_rules([
            ['.', '', 'fasta序列提取结果文件', 0],
        ])
        super(FastaFromIdWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
            "id": "tw_fasta",
            "type":"workflow",
            "name": "tw_fasta_from_id",
            "instant": False,
            "options": {
                "fasta_file": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/test/t.fasta",
                "id_file": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/test/t.ids",
                    }
        }
    wsheet = Sheet(data=data)
    wf = FastaFromIdWorkflow(wsheet)
    wf.run()

