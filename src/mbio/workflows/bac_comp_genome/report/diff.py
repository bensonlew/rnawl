# -*- coding: utf-8 -*- 
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import pandas as pd
import json
import os


class DiffWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(DiffWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'test', 'type': 'string'},  # 检验类型
            {'name': 'testtype', 'type': 'string', 'default': 'two.side'},
            {'name': 'samp1', 'type': 'string', 'default': ''},  # 两样本检验所选样本
            {'name': 'samp2', 'type': 'string', 'default': ''},  # 两样本检验所选样本
            #{'name': 'gfile', 'type': 'infile', 'format': 'meta.otu.group_table'},  # 分组文件
            {'name': 'gfile', 'type': 'string', 'default': '' },  # group_detail字符串
            {'name': 'splist', 'type': 'string', 'default': '' },  #  core pan 分析中pangenome分类比较，选择的样本
            {'name': 'coverage', 'type': 'float', 'default': 0.95},
            {'name': 'norm', 'type': 'string', 'default': 'T'},
            {'name': 'annotable', 'type': 'infile', 'format': 'sample.data_dir'},
            {'name': 'functiontype', 'type': 'string', 'default': ''},
            {'name': 'level', 'type': 'string'},
            {'name': 'corepan', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'pancat', 'type': 'string', 'default': ''},
            {'name': 'selectedcat', 'type': 'string', 'default': 'ALL'},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'main_name', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'status_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

        self.parse_annot = self.add_tool('bac_comp_genome.comp_function')
        self._pasre_out = 'annot_parse.xls'
        self.diff = self.add_tool('bac_comp_genome.diff_statistical')
        self.group_num = None
    
    def check_options(self):
        self.annotable = ''
        for f in os.listdir(self.option('annotable').prop['path']):
            if f.endswith('_anno.xls'):
                self.annotable = os.path.join(self.option('annotable').prop['path'], f)
        if not self.annotable:
            raise OptionError('找不到注释文件 {}'.format(self.option('annotable').path))

    def run(self):
        self.parse_annot.on('end', self.run_diff)
        self.diff.on('end', self.set_db)
        self.set_gfile()
        self.run_parse()
        super(DiffWorkflow, self).run()

    def set_gfile(self):
        self.gfile = self.work_dir + '/gfile.txt'
        self.logger.info(self.gfile)
        with open(self.gfile, 'w') as gf:
            gf.write('#name\tgroup_name\n')
            if self.option('gfile'):
                group_detail = json.loads(self.option('gfile'))
                self.group_num = len(group_detail)
                for k in sorted(group_detail):
                    for n in group_detail[k]:
                        gf.write('{}\t{}\n'.format(n, k))
            elif self.option('splist'):
                [gf.write('{}\t-\n'.format(sp)) for sp in self.option('splist').split(',')]
            else:
                gf.write('{}\t-\n{}\t-\n'.format(self.option('samp1'), self.option('samp2')))

    def run_parse(self):
        options = {
            'annotable': self.annotable,
            'functiontype':  self.option('functiontype'),
            'level': self.option('level'),
            'splist': self.gfile,
            'output': self._pasre_out,
            'abundance': 'T',
            'corepan': self.option('corepan').path,
            'pancat': self.option('pancat'),
            'selectedcat': self.option('selectedcat'),
        }
        self.parse_annot.set_options(options)
        self.parse_annot.run()

    def run_diff(self):
        gfile = self.gfile
        if len(self.option('selectedcat').split(',')) > 1:
            samp1, samp2 = self.option('selectedcat').split(',')
            self.option('samp1', samp1)
            self.option('samp2', samp2)
            gfile = ''
        options = {
            'test': self.option('test'),
            'testtype': self.option('testtype'),
            'intable': self.parse_annot.work_dir + '/' + self._pasre_out,
            'samp1': self.option('samp1'),
            'samp2': self.option('samp2'),
            'gfile': gfile,
            'coverage': self.option('coverage'),
            'norm': self.option('norm')
        }
        self.diff.set_options(options)
        self.diff.run()

    def set_db(self):
        flpath = self.diff.output_dir + '/overall_result.xls'
        try:
            self.link(flpath, 'output/diff_analysis.xls')
        except Exception as e:
            self.logger.info('workflow中，设置结果目录失败{}' % e)
        self.logger.info('开始导表！')

        diff_api = self.api.api('bac_comp_genome.diff_comp')
        #renames = {'corrected_pvalue': 'qvalue', 'lowerCI': 'lowerci', 'upperCI': 'upperci'}
        renames = {'corrected_pvalue': 'qvalue'}
        if self.option('functiontype').upper() == 'KEGG':
            diff_api.update_table(self.option('main_name'), self.option('main_id'),
                                  data={'graph_dir': self.config.SOFTWARE_DIR + '/database/KEGG/map_html/'})
        diff_api.add_detail(flpath, self.option('main_name') + '_detail',
                            self.option('main_id'), self.option('main_name') + '_id',
                            mongo_keys = renames)

        with open(flpath, 'r') as r:
            header = r.readline().strip('\t\r\n').split('\t')
        if self.group_num == 2:
            columns = ['function', 'lowerCI', 'upperCI', 'pvalue'] +\
                    filter(lambda x: x.endswith('mean'), header)
        else:
            if self.option('samp1'):
                columns = ['function', 'lowerCI', 'upperCI', 'pvalue'] + filter(lambda x: x.endswith('propotion'), header)
            else:
                columns = ['function', 'l', 'n', 'pvalue'] + filter(lambda x: x.endswith('mean'), header)
        diff_api.add_detail(flpath, self.option('main_name') + '_plot',
                              self.option('main_id'), self.option('main_name') + '_id',
                              columns, mongo_keys = renames)
        # update sg_status
        status_info = json.loads(self.option('status_info'))
        diff_api.update_sg_status(self.option('main_id'), status_info)
        self.end()

    def end(self):
        # 设置结果目录
        repaths = [
            ['diff_analysis.xls', 'xls', '功能差异分析结果', 0, '']
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)

        super(DiffWorkflow, self).end()
