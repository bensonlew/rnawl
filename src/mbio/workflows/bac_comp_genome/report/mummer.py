# -*- coding: utf-8 -*- 
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
from biocluster.option import OptionError
import json
import os


class MummerWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(MummerWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'mummer', 'type': 'string', 'default': 'nucmer'},  # nucmer promer
            {'name': 'ref', 'type': 'string', 'default': ''},  # 指定的参考基因组
            {'name': 'samples', 'type': 'string'},  # 逗号隔开的样本名
            {'name': 'seq_dir', 'type': 'infile', 'format': 'bac_comp_genome.input_dir'},
            #{'name': 'seq_dir', 'type': 'string',},
            {'name': 'super', 'type': 'bool', 'default': True},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'main_name', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'status_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.mummer = self.add_tool('bac_comp_genome.mummer')

    def run(self):
        self.mummer.on('end', self.set_db)
        self.run_mummer()

        super(MummerWorkflow, self).run()

    def run_mummer(self):
        options = {
            'mummer': self.option('mummer'),
            'ref': self.option('ref'),
            'samples': self.option('samples'),
            'seq_dir': self.option('seq_dir').prop['path'],
            'super': self.option('super')
        }
        self.mummer.set_options(options)
        self.mummer.run()

    def set_db(self):
        mm_api = self.api.api('bac_comp_genome.diff_comp')
        self.logger.info('开始导表！')
        filepath = os.path.join(self.mummer.output_dir, 'mummer.output.xls')
        mm_api.add_detail(filepath, self.option('main_name') + '_detail',
                              self.option('main_id'), self.option('main_name') + '_id'
                              )
        mongo_keys = {0: 'chr', 1: 'len'}
        for sp in os.listdir(self.mummer.work_dir):
            if sp.endswith('.seq_len.xls'):
                s = sp.split('.seq_len.xls')[0]
                fp = os.path.join(self.mummer.work_dir, sp)
                mm_api.add_detail(fp, self.option('main_name') + '_len',
                                  self.option('main_id'), self.option('main_name') + '_id', mongo_keys=mongo_keys,
                                  tag_key=['sp',], tag_value=[s,], header=False)
        self.logger.info('导表结束！')
        self.link(filepath, 'output/syteney.xls')
        self.end()

    def end(self):
        repaths = [
            ['syteney.xls', 'xls', '共性性分析结果', 0, '']        
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)

        super(MummerWorkflow, self).end()
