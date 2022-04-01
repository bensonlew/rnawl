# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
import json
import os


class FunctionEnrichmentWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FunctionEnrichmentWorkflow, self).__init__(wsheet_object)
        options = self.interface()
        self.add_option(options)
        self.set_options(self._sheet.options())

        self.enrich = self.add_tool('bac_comp_genome.function_enrichment') # 修改tool路径
        self.graph = self.add_tool('bac_comp_genome.kegg_graph_info')

    def check_options(self):
        self.annotable = ''
        for f in os.listdir(self.option('annotable').prop['path']):
            if f.endswith('_anno.xls'):
                self.annotable = os.path.join(self.option('annotable').prop['path'], f)
        if not self.annotable:
            raise OptionError('找不到注释文件 {}'.format(self.option('annotable').path))

    def interface(self):
        opts = [
            {'name': 'annotable', 'type': 'infile', 'format': 'sample.data_dir'},
            {'name': 'samples', 'type': 'string'},
            {'name': 'testset', 'type': 'string'},
            {'name': 'corrected', 'type': 'string'},
            {'name': 'functype', 'type': 'string'},
            {'name': 'levelenriched', 'type': 'string'},
            {'name': 'graph', 'type': 'bool', 'default': False},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string', 'default': ''},
            {'name': 'main_name', 'type': 'string', 'default': ''},
        ]
        return opts

    def run_first(self):
        '''
        运行没前置依赖的步骤
        '''
        self.run_enrich()

    def run(self):
        rls = [self.enrich,]
        if self.option('graph') and self.option('functype').upper() == 'KEGG':
            rls.append(self.graph)
            self.get_graph()
        self.on_rely(rls, self.set_db)

        self.run_first()
        super(FunctionEnrichmentWorkflow, self).run()

    def run_enrich(self):
        options = {
            'annotable': self.annotable, 
            'samples': self.option('samples'),
            'testset': self.option('testset'),
            'enrichout': self.option('functype') + '_enrichment.xls',
            'corrected': self.option('corrected'),
            'levelenriched': self.option('levelenriched'),
            'formated': 'N',
            'functype': self.option('functype')
        }
        self.enrich.set_options(options)

        self.enrich.run()

    def get_graph(self):
        options = {
            'annotable': self.annotable,
        }
        self.graph.set_options(options)

        self.graph.run()

    def set_db(self):
        cm_api = self.api.api('bac_comp_genome.diff_comp')

        self.logger.info('开始导表')
        if self.option('functype').upper() == 'KEGG':
            cm_api.update_table(self.option('main_name'), self.option('main_id'),
                                  data={'graph_dir': self.config.SOFTWARE_DIR + '/database/KEGG/map_html/'})

        enrich_out = self.enrich.work_dir + '/' + self.option('functype') + '_enrichment.xls'
        files = [[enrich_out, self.option('main_name') + '_detail'],]

        if self.option('graph') and self.option('functype').upper() == 'KEGG':
            graph_out = self.graph.work_dir + '/' + 'kegg_graph_info.xls'

            files.append([graph_out, self.option('main_name') + '_graph'])

        for f in files:
            with open(f[0], 'r') as r:
                header = r.readline().replace('.', '_')
            mongo_key = header.strip('\n\t').split('\t')
            cm_api.add_detail(f[0], f[1], self.option('main_id'), 
                              self.option('main_name') + '_id', header_lower=True)
            try:
                self.link(f[0])
            except Exception as e:
                self.logger.info('workflow中，设置结果目录失败{}' % e)
        self.end()

    def end(self):
        regexps = [
            [r'.*_enrichment.xls', 'xls', self.option('functype') + '富集分析结果', 0, ''],
        ]
        repaths = []
        if self.option('functype').upper() == 'KEGG':
            repaths = [
                ['kegg_graph_info.xls', 'xls', 'kegg HTML图位置信息', 0, ''],        
            ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

        super(FunctionEnrichmentWorkflow, self).end()
