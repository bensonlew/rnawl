# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
from biocluster.tool import Tool
import pandas as pd
import unittest

class ExpVennAgent(Agent):
    '''
    require express_matrix (file), group_table (file) and threshold (float)
    '''
    def __init__(self, parent):
        super(ExpVennAgent, self).__init__(parent)
        options = [
            {'name': 'express_matrix', 'type': 'infile', 'format': 'small_rna.express_matrix'},
            {'name': 'group_table', 'type': 'infile', 'format': 'small_rna.group_table'},
            {'name': 'threshold', 'type': 'float', 'default': 1.0}
        ]
        self.add_option(options)

    def check_options(self):
        try:
            threshold = float(self.option('threshold'))
            self.logger.debug('type of incoming threshold is {}, pass'.format(type(threshold)))
        except:
            self.logger.debug('type of incoming threshold is {}, abord'.format(type(threshold)))
            raise OptionError('type of threshold can not be converted to float')
        if threshold <= 0:
            self.logger.debug('value of incoming threshold is {}, abord'.format(threshold))
            raise OptionError('value of threshold must be greater than zero')

    def set_resource(self):
        self._cpu = 1
        file_size = os.path.getsize(self.option('express_matrix').prop['path'])
        self._memory = '{}G'.format(int(float(file_size)/1024**3) + 4)

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     ['.', '', 'exp_venn_tool_output_dir'],
        # ])
        super(ExpVennAgent, self).end()

class ExpVennTool(Tool):
    '''
    obtain exp_venn result, write data to venn_graph.xls
    '''
    def __init__(self, config):
        super(ExpVennTool, self).__init__(config)
        self.venn_graph = 'venn_graph.xls'

    def express_venn(self):
        self.logger.info('start express_venn at tool')
        threshold = float(self.option('threshold'))
        group_dict = self.option('group_table').prop['group_dict']
        exp_table = pd.read_table(self.option('express_matrix').prop['path'], index_col=0, header=0)
        all_sets = dict()
        self.logger.debug('write venn_graph to {}'.format(os.path.join(self.work_dir, self.venn_graph)))
        f = open(os.path.join(self.work_dir, self.venn_graph), 'w')
        f.write('#set_name\tset_detail\n')
        for each_group in group_dict:
            samples = group_dict[each_group]
            expressed = exp_table.index[exp_table.loc[:, samples].mean(axis=1) >= threshold]
            all_sets[each_group] = list(expressed)
            f.write(each_group+'\t'+','.join(expressed)+'\n')
        f.close()
        self.logger.info('finish express_venn at tool')

    def set_output(self):
        self.logger.info('start set_output at tool')
        source = os.path.join(self.work_dir, self.venn_graph)
        link_name = os.path.join(self.output_dir, self.venn_graph)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at tool')

    def run(self):
        super(ExpVennTool, self).run()
        self.express_venn()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):

        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import random

        exp_matrix_list = ['/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_count.xls',
                           '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls',
                           '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_count.xls',
                           '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_norm.xls']

        for exp_matrix in exp_matrix_list:

            data = {
                'id': 'exp_venn_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'tool',
                'name': 'small_rna.exp_venn',
                'instant': True,
                'options': {
                    'express_matrix': exp_matrix,
                    'group_table': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/group.txt',
                    'threshold': 1.0,
                }
            }

            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

if __name__ == '__main__':
    unittest.main()