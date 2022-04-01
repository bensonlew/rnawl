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


class Fasta2sqnWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(Fasta2sqnWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'temp_file', 'type': 'infile', 'format': 'tool_lab.no_empty'},
            {'name': 'seq_file', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'structure_type', 'type': 'string'},
            {'name': 'structure_file', 'type': 'infile', 'format': 'tool_lab.no_empty'},
            {'name': 'annot_file', 'type': 'infile', 'format':'sequence.profile_table'},
            {'name': 'species', 'type': 'string'},
            {'name': 'lab_name', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.f2s = self.add_tool('tool_lab.fasta2sqn')

    def run(self):
        self.f2s.on('end', self.end)
        self.run_f2s()
        super(Fasta2sqnWorkflow, self).run()

    def run_f2s(self):
        opts = {
            'temp_file': self.option('temp_file').path,
            'seq_file': self.option('seq_file').path,
            'structure_type': self.option('structure_type'),
            'structure_file': self.option('structure_file').path,
            'annot_file': self.option('annot_file').path,
            'species': self.option('species'),
            'lab_name': self.option('lab_name'),
        }
        if self.option('annot_file').is_set:
            opts['annot_file'] = self.option('annot_file').path
        if self.option('species'):
            opts['species'] = self.option('species')
        if self.option('lab_name'):
            opts['lab_name'] = self.option('lab_name')
        self.f2s.set_options(opts)
        self.f2s.run()

    def end(self):
        self.logger.info('开始上传结果文件')
        sqn_tag = ''
        for f in os.listdir(self.f2s.output_dir):
            if f.endswith('.sqn'):
                self.link(os.path.join(self.f2s.output_dir, f))
                sqn_tag = '1'
        if not sqn_tag:
            self.set_error('运行失败未生成sqn文件')
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', "", "fasta序列转sqn文件夹"],    
        ])
        result_dir.add_regexp_rules([
            [r'.*\.sqn', "sqn", "fasta转sqn结果文件"],
        ])
        self.logger.info('成功上传结果文件')
        super(Fasta2sqnWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        "id": "tt_sqn",
        "type":"workflow",
        "client": "client03",
        "rerun": True,
        "options": {
            "temp_file": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/test/template.sbt",
            "seq_file": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/test/sequence.fsa",
            "structure_type": "tbl",
            "structure_file": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/test/test.tbl",
            }
    }

    wsheet = Sheet(data=data)
    wf = Fasta2sqnWorkflow(wsheet)
    wf.run()

