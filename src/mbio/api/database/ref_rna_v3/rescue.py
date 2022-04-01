# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test_all_exp(self):
        import random
        from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rescue_all_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v2.refrna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = RefrnaTestApiWorkflow(wheet)
        wf.sheet.id = 'i-sanger_132345'
        wf.sheet.project_sn = '12923_5bd272b137663'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False

        output_dir = '/mnt/lustre/users/sanger/sg-users/rescue/i-sanger_132345/Express/ExpAnnalysis'

        api = wf.api.api('ref_rna_v2.all_exp')
        exp_type = 'TPM'
        exp_matrix = {'T': os.path.join(output_dir, 'transcript.{}.matrix.xls'.format(exp_type.lower())),
                      'G': os.path.join(output_dir, 'gene.{}.matrix.xls'.format(exp_type.lower()))}
        quant_method = 'RSEM'
        lib_type = None
        group_dict = {'IBU': ['Ep_IBU1', 'Ep_IBU2', 'Ep_IBU3'], 'MJ': ['Ep_MJ1', 'Ep_MJ2', 'Ep_MJ3']}
        group_id = '5be5c6f028fb4f214bd538ed'
        project_sn = '12923_5bd272b137663'
        task_id = 'i-sanger_132345'
        params = json.dumps({
            'task_id': task_id,
            'submit_location': 'exp_detail',
            'task_type': 2,
            'method': quant_method,
            'exp_type': exp_type
        }, sort_keys=True, separators=(',', ':'))
        exp_ids = dict()
        exp_ids['T'] = api.add_exp(exp_matrix=exp_matrix['T'], quant_method=quant_method, exp_level='T',
                                   lib_type=lib_type, group_dict=group_dict, group_id=group_id,
                                   exp_type=exp_type, add_distribution=False, project_sn=project_sn,
                                   task_id=task_id, params=params)
        exp_ids['G'] = api.add_exp(exp_matrix=exp_matrix['G'], quant_method=quant_method, exp_level='G',
                                   lib_type=lib_type, group_dict=group_dict, group_id=group_id,
                                   exp_type=exp_type, add_distribution=False, project_sn=project_sn,
                                   task_id=task_id, params=params)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_all_exp')])
    unittest.TextTestRunner(verbosity=2).run(suite)
