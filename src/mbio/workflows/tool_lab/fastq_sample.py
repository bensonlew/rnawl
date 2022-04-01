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


class FastqSampleWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(FastqSampleWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'read1', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'read2', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'read', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'reads_type', 'type': 'string', },
            {'name': 'extract_type', 'type': 'string', },
            {'name': 'size_value', 'type': 'int',},
            {'name': 'base_value', 'type': 'int',},
            {'name': 'seq_num_value', 'type': 'int',},
            {'name': 'threshold', 'type': 'float', },
            {'name': 'nthread', 'type': 'int', 'default': 2},
            {'name': 'main_table', 'type': 'string', 'default': 'fastq_sample'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.fq_sp = self.add_tool('tool_lab.fastq_sample')

    def run(self):
        self.fq_sp.on('end', self.set_db)
        self.run_fq_sp()
        super(FastqSampleWorkflow, self).run()

    def run_fq_sp(self):
        vl = self.option('extract_type') + '_value'
        opts = {
            'extract_type': self.option('extract_type'),
            'threshold': self.option(vl),
            'nthread': self.option('nthread'),
        }
        if self.option('read').is_set:
            opts['read'] = self.option('read').prop['path']
        else:
            opts['read1'] = self.option('read1').prop['path']
            opts['read2'] = self.option('read2').prop['path']
        self.fq_sp.set_options(opts)
        self.fq_sp.run()

    def set_db(self):
        #remote_dir = self._sheet.output
        #file_path = os.path.join(remote_dir, self.option('output'))
        #api = self.api.api('tool_lab.api_base')
        #query_dict = {'_id': self.option('main_id')}
        #update_dict = {'result_path': file_path}
        #api.update_db_record(self.option('main_table'), query_dict, update_dict)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.fq_sp.output_dir)
        result_dir.add_regexp_rules([
            [r'.*.fq', 'fq', '提取结果文件', 0],
            [r'.*.fastq', 'fastq', '提取结果文件', 0],
        ])
        result_dir.add_relpath_rules([
            ['.', '', '提取fq文件夹', 0],
        ])
        super(FastqSampleWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet

    data = {
            'id': 'tw_fastq',
            'name': 'test_fastq_sample',
            'type': 'workflow',
            'options': {
                "read1": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/r1.fq",
                "read2": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/r2.fq",
                "extract_type": "seq_num",
                "seq_num_value": 4000,
                }
            }
    wsheet = Sheet(data=data)
    wf = FastqSampleWorkflow(wsheet)
    wf.run()

