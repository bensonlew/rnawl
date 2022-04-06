# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.agent import Agent
import os
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import pandas as pd
from mbio.packages.small_rna.go_graph import draw_GO
import json
import datetime
import glob
import unittest

class GoDagAgent(Agent):
    def __init__(self, parent):
        super(GoDagAgent, self).__init__(parent)
        options = [
            {"name": "go_enrich_detail", "type": "infile", 'format': 'small_rna.common'},
            {"name": "significant_diff", "type": "string", 'default': 'pvalue'},
            {"name": "significant_value", "type": "string", 'default': '0.05'},
            {"name": "top_num", "type": "int", 'default': 10},
            {"name": "go_list", "type": "string", "default": None},
        ]
        self.add_option(options)
        self.step.add_steps('go_dag')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.go_dag.start()
        self.step.update()

    def stepfinish(self):
        self.step.go_dag.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self._options.keys():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        if not os.path.exists(self.option('go_enrich_detail').prop['path']):
            raise OptionError('{} not exist'.format(self.option('go_enrich_detail').prop['path']))
        if self.option('go_list') is None and self.option('significant_diff') is None:
            raise OptionError('find both None in significant_diff and go_list')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        super(GoDagAgent,self).end()

class GoDagTool(Tool):
    def __init__(self,config):
        super(GoDagTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.go_graph = os.path.join(self.config.PACKAGE_DIR, 'small_rna/go_graph.py')
        self.obo = os.path.join(self.config.SOFTWARE_DIR, 'database/GO/go-basic.obo')

    def generate_go_dict(self):
        go_p_dict = dict()
        if self.option('significant_diff').lower() == 'pvalue':
            by_fdr = False
        else:
            by_fdr = True
        go_enrich_geneset_sort = self.sort_p(self.option('go_enrich_detail').prop['path'], by_fdr)
        with open(go_enrich_geneset_sort) as f:
            _ = f.readline()
            if self.option('go_list') is None:
                self.logger.debug('')
                count = 0
                for line in f:
                    if self.option('significant_diff').lower() == 'pvalue':
                        significant_diff = line.split('\t')[6]
                    else:
                        significant_diff = line.split('\t')[9]
                    if count < self.option('top_num'):
                        if float(significant_diff) < float(self.option('significant_value')):
                            go_p_dict[line.split('\t')[0]] = float(significant_diff)
                    else:
                        break
                    count = count + 1
            else:
                self.logger.debug('')
                for line in f:
                    for items in self.option('go_list').split(';'):
                        if line.startswith(items):
                            go_p_dict[line.split('\t')[0]] = float(line.split('\t')[9])
        return go_p_dict

    @staticmethod
    def sort_p(tabular_file, by_fdr=False):
        df = pd.read_table(tabular_file)
        if by_fdr:
            df.sort_values(by='p_fdr_bh', inplace=True)
        else:
            df.sort_values(by='p_uncorrected', inplace=True)
        sort_file = '{}.sort'.format(tabular_file)
        df.to_csv(sort_file, sep='\t', index=None)
        return sort_file

    def run_dag(self):
        go_p_dict = self.generate_go_dict()
        if len(go_p_dict) == 0:
            if self.option('significant_diff').lower() == 'pvalue':
                self.set_error('no data pass the given threshold, abord')
            else:
                self.set_error('no data pass the given threshold, please choose pvalue as filter type')
        draw_GO(go_p_dict, out='go_dag', obo=self.obo)

    def set_output(self):
        go_dag = glob.glob(self.work_dir + '/go_dag.*')
        for each in go_dag:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(GoDagTool, self).run()
        self.run_dag()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_default(self):

        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'go_dag_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'small_rna.geneset.go_dag',
            'instant': False,
            'options': {
                'go_enrich_detail': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/geneset_go_dag/go_enrich_geneset.xls',
                'significant_diff': 'pvalue',
                'significant_value': '0.05',
                'top_num': 10,
            }
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_go_list(self):

        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'go_dag_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'small_rna.geneset.go_dag',
            'instant': False,
            'options': {
                'go_enrich_detail': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/geneset_go_dag/go_enrich_geneset.xls',
                'go_list': 'GO:0044464;GO:0044424;GO:0005488;GO:0009987;GO:0065007;GO:0043226;GO:0044699;GO:0043227;GO:0043229;GO:0044763;GO:0043231;GO:0044444;GO:0044422;GO:0044446;GO:0005515;GO:0019222;GO:0043167;GO:0031323;GO:0060255;GO:0080090',
            }
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()