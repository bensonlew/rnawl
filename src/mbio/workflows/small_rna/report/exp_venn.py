# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.workflow import Workflow
import os
import time
import unittest
import json

class ExpVennWorkflow(Workflow):
    '''
    transfer exp_matrix (file:str), group (file:str) and threshold (float:str) to tool
    obtain venn_graph.xls from tool
    deliver venn_graph.xls to api according to main_id of sg_exp_venn
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpVennWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'string'},
            {'name': 'group_dict', 'type': 'string'},  # to_file require, no use here
            {'name': 'seq_type', 'type': 'string'},  # to_file require, no use here
            {'name': 'group', 'type': 'string'},
            {'name': 'threshold', 'type': 'string', 'default': '1.0'},
            {'name': 'venn_main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool('small_rna.exp_venn')
        self.add_exp = self.api.api('small_rna.all_exp')

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self.sheet.options():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''
        self.tool.on('end', self.set_output)
        self.run_tool()
        super(ExpVennWorkflow, self).run()

    def run_tool(self):
        '''
        set options for tool
        set trigger point
        '''
        options = {
            'express_matrix': self.option('exp_matrix'),
            'group_table': self.option('group'),
            'threshold': self.option('threshold'),
        }
        self.tool.set_options(options)
        self.tool.run()

    def set_output(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = os.path.join(self.tool.output_dir, 'venn_graph.xls')
        link_name = os.path.join(self.output_dir, 'venn_graph.xls')
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        '''
        export result in output_dir of workflow to api
        '''
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        venn_graph = os.path.join(self.output_dir, 'venn_graph.xls')
        self.add_exp.add_exp_venn(venn_graph, main_id=self.option('venn_main_id'))
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        '''
        declare output_dir and call super class method to end workflow
        '''
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     ['.', '', '样本间Venn分析结果文件', 0],
        #     ['venn_graph.xls', '', 'Venn分析对应详情表', 0],
        # ])
        super(ExpVennWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_add_known_tpm(self):

        from mbio.workflows.small_rna.report.exp_venn import ExpVennWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.report.exp_venn',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/all_miR_norm.xls',
                'seq_type': 'all',
                'group': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/group.txt',
                'group_dict': json.dumps({'GH': ['GH1', 'GH2'], 'NFAs': ['NFAs1', 'NFAs2'], 'Normal': ['Normal1', 'Normal2'], 'PRL': ['PRL1', 'PRL2']}),
                'threshold': '1.0',
                'venn_main_id': None
            }
        }

        wsheet = Sheet(data=data)
        wf = ExpVennWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()