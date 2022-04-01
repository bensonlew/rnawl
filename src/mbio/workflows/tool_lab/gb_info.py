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


class GbInfoWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(GbInfoWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'gb_file', 'type': 'infile', 'format': 'gene_structure.gbk'},
            {'name': 'extract', 'type': 'string', 'default': ''},
            {'name': 'cds', 'type': 'int', 'default': 0},
            {'name': 'genome', 'type': 'int', 'default': 0},
            {'name': 'prot', 'type': 'int', 'default': 0},
            {'name': 'gff', 'type': 'int', 'default': 0},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.gb_info = self.add_tool('tool_lab.gb_info')

    def run(self):
        self.gb_info.on('end', self.end)
        self.run_gb_info()
        super(GbInfoWorkflow, self).run()

    def run_gb_info(self):
        selected = self.option('extract').split(',')
        [self.option(k, 1) for k in selected]
        opts = {
            'gb_file': self.option('gb_file').path,
            'cds': self.option('cds'),
            'prot': self.option('prot'),
            'gff': self.option('gff'),
            'genome': self.option('genome')
        }
        self.gb_info.set_options(opts)
        self.gb_info.run()
    
    def end(self):
        self.logger.info('开始上传结果文件')
        relpath = [
            ['.', "", "GenBank文件信息提取结果文件夹"],
        ]
        regpath = [
            [r'.*fasta', "fasta", "基因组文件"],
            [r'.*cds', "cds", "cds序列"],
            [r'.*faa', "faa", "faa序列"],
            [r'.*gff', "gff", "gff文件"],
        ]
        up_dir = self.add_upload_dir(self.gb_info.output_dir)
        up_dir.add_relpath_rules(relpath)
        up_dir.add_regexp_rules(regpath)
        self.logger.info('成功上传结果文件')
        super(GbInfoWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        "id": "tt_gb",
        "type": "workflow",
        "name": "tw_gb_info",
        "instant": False,
        "options": {
                "gb_file": "/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/coding/package/EMBL_GBK2gff/input/NC_002695.gbk",
                "extract": "cds,genome",
                #"cds": 1,
                #"prot": 1,
                #"gff": 1,
                #"genome": 1,
            }
    }
    wsheet = Sheet(data=data)
    wf = GbInfoWorkflow(wsheet)
    wf.run()
