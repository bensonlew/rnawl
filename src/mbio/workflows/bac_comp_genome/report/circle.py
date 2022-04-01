# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
import json
import os


class CircleWorkflow(Workflow):
    def __init__(self, wsheet, **kwargs):
        super(CircleWorkflow, self).__init__(wsheet, **kwargs)
        options = [
            {'name': 'ref', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},
            {'name': 'seq_dir', 'type': 'infile', 'format': 'bac_comp_genome.input_dir'},
            {'name': 'window', 'type': 'int', 'default': 10000},
            {'name': 'step', 'type': 'int', 'default': 10000},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'main_name', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'status_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(wsheet.options())
        self.gc_infos = self.add_tool('bac_comp_genome.gc_infos')
        self.mummer = self.add_tool('bac_comp_genome.mummer')

    def run(self):
        self.on_rely([self.gc_infos, self.mummer], self.set_db)
        self.run_gc_infos()
        self.run_mummer()
        super(CircleWorkflow, self).run()

    def run_gc_infos(self):
        options = {
            'ref': self.option('ref'),
            'seq_dir': self.option('seq_dir'),
            'window': self.option('window'),
            'step': self.option('step')
        }
        self.gc_infos.set_options(options)
        self.gc_infos.run()

    def run_mummer(self):
        options = {
            'ref': self.option('ref'),
            'samples': self.option('samples'),
            'seq_dir': self.option('seq_dir'),
            'circle_mode': True
        }
        self.mummer.set_options(options)
        self.mummer.run()

    def set_db(self):
        my_api = self.api.api('bac_comp_genome.diff_comp')
        self.logger.info('开始导表')
        prefix = self.option('ref')
        gc_content = os.path.join(self.gc_infos.work_dir, prefix + '.gc_content.xls')
        gc_skew = os.path.join(self.gc_infos.work_dir, prefix + '.gc_skew.xls')
        gc_mongo = {0: 'location', 1: 'start', 2: 'end', 3: 'value'}
        my_api.add_detail(gc_content, self.option('main_name') + '_detail',
                          self.option('main_id'), self.option('main_name') + '_id', mongo_keys=gc_mongo,
                          tag_key=['type', 'sample'], tag_value=['gc_content', self.option('ref')],
                          header=False)
        my_api.add_detail(gc_skew, self.option('main_name') + '_detail',
                          self.option('main_id'), self.option('main_name') + '_id', mongo_keys=gc_mongo,
                          tag_key=['type', 'sample'], tag_value=['gc_skew', self.option('ref')],
                          header=False)
        
        mm_file = os.path.join(self.mummer.output_dir, 'mummer.output.xls')
        col = ['query', 'ref_scaf', 'ref_start', 'ref_end', 'identity']
        rename = ['sample', 'location', 'start', 'end', 'value']
        mongo_keys = dict(zip(col, rename))
        my_api.add_detail(mm_file, self.option('main_name') + '_detail',
                          self.option('main_id'), self.option('main_name') + '_id', columns=col,
                          mongo_keys=mongo_keys, tag_key=['type', ], tag_value=['query', ])
                    
        self.link(gc_skew, 'output/gc_skew.xls')
        self.link(gc_content, 'output/gc_content.xls')
        self.link(mm_file, 'output/circle.xls')
        self.end()

    def end(self):
        repaths = [
            ['gc_skew.xls', 'xls', 'gc skew数据表', 0, ''],
            ['gc_content.xls', 'xls', 'gc content数据表', 0, ''],
            ['circle.xls', 'xls', '参考基因组和其它基因的同源区域数据表', 0, '']
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)

        super(CircleWorkflow, self).end()
