# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import glob
import pandas as pd
import unittest

class SmallRnaTestExpWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SmallRnaTestExpWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'known_matrix', 'type': 'string', 'default': ''},
            {'name': 'novel_matrix', 'type': 'string', 'default': ''},
            {'name': 'project_sn', 'type': 'string', 'default': 'small_rna'},
            {'name': 'task_id', 'type': 'string', 'default': 'small_rna'},
            {'name': 'exp_type', 'type': 'string', 'default': 'tpm'},
            {'name': 'group', 'type': 'infile', 'format': 'small_rna.group_table'},
            {'name': 'group_id', 'type': 'string', 'default': ''},
            {'name': 'update_info', 'type': 'string', 'default': ''}
        ]
        self.task_id = self.sheet.id
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.all_exp = self.api.api('small_rna.all_exp')

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for opt in self.sheet.options():
            self.logger.debug('{} - {}'.format(opt, self.option(opt)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        self.start_listener() # start_listener() is essential for end()
        self.fire('start')
        self.set_output()
        self.end() # tigger end() exactly which is depend on start_listener()

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        # make exp_matrix
        all_exp_pd_lst = list()
        for each in [self.option('known_matrix'), self.option('novel_matrix')]:
            try:
                all_exp_pd_lst.append(self.all_exp.process_exp_matrix(each))
            except Exception as e:
                self.set_error('{} can not be processed by pd.DataFrame with exception: {}'.format(each, e))
        all_exp_merged_pd = pd.concat(all_exp_pd_lst, axis=0)
        if not self.option('exp_type'):
            self.set_error('option exp_type is {}, abord'.format(self.option('exp_type')))
        elif self.option('exp_type') == 'tpm':
            self.exp_matrix = os.path.join(self.work_dir, 'all_miR_tpm.xls')
        elif self.option('exp_type') == 'count':
            self.exp_matrix = os.path.join(self.work_dir, 'all_miR_count.xls')
        else:
            self.set_error('option exp_type is {}, abord'.format(self.option('exp_type')))
        all_exp_merged_pd.to_csv(self.exp_matrix, sep='\t', header=True, index=True)
        # link file
        for source in glob.glob(os.path.join(self.work_dir, '*.xls')):
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        known_matrix = self.option('known_matrix')
        novel_matrix = self.option('novel_matrix')
        project_sn = self.option('project_sn')
        task_id = self.option('task_id')
        exp_type = self.option('exp_type')
        group_dict = self.option('group').prop['group_dict']
        group_id = self.option('group_id')
        params = {
            'task_id': task_id,
            'submit_location': 'workflow',
            'task_type': 2,
            'group_id': group_id
        }
        main_id = self.all_exp.add_exp(exp_matrix=known_matrix, project_sn=project_sn, task_id=task_id,
                                       exp_type=exp_type, is_novel=False, params=params)
        exp_id = self.all_exp.add_exp(exp_matrix=novel_matrix, project_sn=project_sn, task_id=task_id,
                                       exp_type=exp_type, is_novel=True, params=params, main_id=main_id)
        if exp_type == 'tpm':
            self.logger.info('start set_db for exp distribution while exp_type is tpm')
            self.all_exp.add_distribution(exp_matrix=self.exp_matrix, group_dict=group_dict,
                                          project_sn=project_sn, task_id=task_id,
                                          exp_id=exp_id, seq_type='all', exp_type='tpm', params=params)
            self.logger.info('finish set_db for exp distribution while exp_type is tpm')
        else:
            self.logger.info('skip set_db at for exp distribution while exp_type is count')
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))

    def end(self):
        super(SmallRnaTestExpWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_tpm(self):

        from mbio.workflows.small_rna.small_rna_test_exp import SmallRnaTestExpWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.small_rna_test_exp',
            'options': {
                'known_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls',
                'novel_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_norm.xls',
                'project_sn': 'small_rna',
                'task_id': 'small_rna',
                'exp_type': 'tpm',
                'group': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/group.txt',
                'group_id': '5be24abac82c5131355dcd00',
            }
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestExpWorkflow(wsheet)
        wf.run()

    def test_count(self):

        from mbio.workflows.small_rna.small_rna_test_exp import SmallRnaTestExpWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.small_rna_test_exp',
            'options': {
                'known_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_count.xls',
                'novel_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_count.xls',
                'project_sn': 'small_rna',
                'task_id': 'small_rna',
                'exp_type': 'count',
                'group': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/group.txt',
                'group_id': '5be24abac82c5131355dcd00',
            }
        }
        wsheet = Sheet(data=data)
        wf = SmallRnaTestExpWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()