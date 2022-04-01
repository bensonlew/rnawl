# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
import json
import os


class GeneclusterCompareWorkflow(Workflow):
    def __init__(self, wsheet, **kwargs):
        super(GeneclusterCompareWorkflow, self).__init__(wsheet, **kwargs)
        options = [
            {'name': 'gene_dir', 'type': 'infile', 'format', 'bac_comp_genome.input_dir'},
            {'name': 'ref', 'type': 'string'},
            {'name': 'seq_dir', 'type': 'infile', 'format', 'bac_comp_genome.input_dir'},
            {'name': 'region', 'type': 'string'},
            {'name': 'cu_names', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'main_name', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'status_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(wsheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        # self.runs字典存放每个簇名对应的module对象
        self.runs = {}

    def run(self):
        ref = json.loads(self.option('ref'))
        region = json.loads(self.option('region'))
        samples = json.loads(self.option('samples'))
        cu_names = json.loads(self.option('cu_names'))
        self.logger.info('### samples')
        self.logger.info(samples)
        # 按簇名add_module
        for cu in cu_names:
            m = self.add_module('bac_comp_genome.genecluster_compare')
            self.runs[cu] = m
            # 保证参考基因组在samples中的第一个
            sps = set(samples[cu].split(';'))
            sps.remove(ref[cu])
            sps = [ref[cu], ] + list(sps)
            samples[cu] = ','.join(sps)
            region[cu] = ','.join(region[cu].split('--'))
        
        self.on_rely(self.runs.values(), self.set_db)
        for cu in cu_names:
            self.one_cu_run(self.runs[cu], region[cu], cu, samples[cu])

        super(GeneclusterCompareWorkflow, self).run()

    def one_cu_run(self, m, region, cu_name, samples):
        opts = {
            'gene_dir': self.option('gene_dir'),
            'seq_dir': self.option('seq_dir'),
            'region': region,
            'cu_name': cu_name,
            'samples': samples,
        }
        m.set_options(opts)
        m.run()

    def set_db(self):
        my_api = self.api.api('bac_comp_genome.diff_comp')
        self.logger.info('开始导表')
        for cu, md in self.runs.items():
            flpath = os.path.join(md.work_dir, cu + '_genecluster_compare.xls')
            self.link(flpath)
            my_api.add_detail(flpath, self.option('main_name') + '_detail',
                              self.option('main_id'), self.option('main_name') + '_id',
                              tag_key = ['cu_name',], tag_value=[cu,]
                            )
 
        self.end()

    def end(self):
        regexps = [
            [r'*_genecluster_compare.xls', 'xls', '基因组比较结果表', 0, '']  
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_regexp_rules(regexps)

        super(GeneclusterCompareWorkflow, self).end()
