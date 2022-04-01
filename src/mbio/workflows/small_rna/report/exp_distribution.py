# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.workflow import Workflow
import os
import time
import unittest
import json

class ExpDistributionWorkflow(Workflow):
    '''
    TODO: description
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpDistributionWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'string'},
            {'name': 'group_dict', 'type': 'string'},  # to_file require, no use here
            {'name': 'seq_type', 'type': 'string'},  # to_file require, no use here
            {'name': 'graph_main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
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
        self.start_listener()
        self.fire('start')
        self.set_db()
        self.end()

    def set_db(self):
        '''
        export result in output_dir of workflow to api
        '''
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        self.add_exp.add_distribution(exp_matrix=self.option('exp_matrix'),
                                      group_dict=json.loads(self.option('group_dict')),
                                      main_id=self.option('graph_main_id'))
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))

    def end(self):
        '''
        declare output_dir and call super class method to end workflow
        '''
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'exp_distribution_workflow_output_dir'],
        ])
        super(ExpDistributionWorkflow, self).end()