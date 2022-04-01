# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
import json
import os


class CompWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(CompWorkflow, self).__init__(wsheet)
        options = [
                {'name': 'annotable', 'type': 'infile', 'format': 'sample.data_dir'},  # 注释文件
            {'name': 'functiontype', 'type': 'string', 'default': ''},  # 注释类型
            {'name': 'level', 'type': 'string', 'default': ''},  # 选择注释水平
            #{'name': 'gfile', 'type': 'infile', 'format': 'meta.otu.group_table'},  # 样本分组
            {'name': 'gfile', 'type': 'string', 'default': ''},  # 样本分组
            {'name': 'groups', 'type': 'string'},  # 选择的分组
            {'name': 'corepan', 'type': 'infile', 'format': 'sequence.profile_table'},  # pangenome结果
            {'name': 'selectedcat', 'type': 'string', 'default': 'ALL'},
            {'name': 'pancat', 'type': 'string', 'default': ''},  # pan分类方案
            {'name': 'result_type', 'type': 'string', 'default': ''},  # 组内合并选项
            {'name': 'main_id', 'type': 'string', 'default': ''}, # 比较分析主表ID
            {'name': 'main_name', 'type': 'string', 'default': ''}, # 比较分析详情表
            {'name': 'update_info', 'type': 'string'},
            {'name': 'status_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.comp = self.add_tool('bac_comp_genome.comp_function') 
        self.outfile = '_function_comparison.xls'

    def check_options(self):
        self.annotable = ''
        for f in os.listdir(self.option('annotable').prop['path']):
            if f.endswith('_anno.xls'):
                self.annotable = os.path.join(self.option('annotable').prop['path'], f)
        if not self.annotable:
            raise OptionError('找不到注释文件 {}'.format(self.option('annotable').path))

    def run(self):
        self.comp.on('end', self.set_db)
        self.set_gfile()
        self.run_comp()
        super(CompWorkflow, self).run()

    def set_gfile(self):
        self.gfile = self.work_dir + '/gfile.txt'
        with open(self.gfile, 'w') as gf:
            gf.write('#name\tgroup_name\n')
            group_detail = json.loads(self.option('gfile'))
            for k in sorted(group_detail):
                for n in group_detail[k]:
                    gf.write('{}\t{}\n'.format(n, k))

    def run_comp(self):
        self.outfile = self.option('functiontype') + self.outfile
        options = {
            'annotable': self.annotable,
            'functiontype': self.option('functiontype'),
            'level': self.option('level'),
            'splist': self.gfile,
            'groups': self.option('groups'),
            'corepan': self.option('corepan').path,
            'pancat': self.option('pancat'),
            'selectedcat': self.option('selectedcat'), 
            'output': self.outfile,
            'result_type': self.option('result_type'),
        }
        self.comp.set_options(options)
        self.comp.run()

    def set_db(self):
        # 设置结果目录
        flpath = self.comp.work_dir + '/' + self.outfile
        try:
            self.link(flpath)
        except Exception as e:
            self.logger.info('workflow中，设置结果目录失败{}' % e)

        self.logger.info('开始导表！')
        cm_api = self.api.api('bac_comp_genome.diff_comp')
        if self.option('functiontype').upper() == 'KEGG':
            cm_api.update_table(self.option('main_name'), self.option('main_id'),
                                  data={'graph_dir': self.config.SOFTWARE_DIR + '/database/KEGG/map_html/'})
        cm_api.add_detail(flpath, self.option('main_name') + '_detail',
                               self.option('main_id'),
                               self.option('main_name') + '_id'
                               )
        status_info = json.loads(self.option('status_info'))
        cm_api.update_sg_status(self.option('main_id'), status_info)
        self.logger.info('导表结束！')

        self.end()

    def end(self):
        repaths = [
            [self.outfile, 'xls', self.option('functiontype') + '统计比较结果表', 0, '']
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)

        super(CompWorkflow, self).end()

